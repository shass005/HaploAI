import os
import glob
import subprocess
import pandas as pd
import numpy as np
import json
import streamlit as st
from streamlit_option_menu import option_menu
import streamlit.components.v1 as components
from scipy.optimize import nnls as scipy_nnls

# Class for processing the target file
class TargetProcessor:
    # initial the paths and create output directories
    def __init__(
        self,
        main_path="plink/ref_pca",
        plink_exe="plink/plink1.9.exe",
        plink2_exe="plink/plink2.exe",
        output_dir="dtc_pcs",
        normalised_dir="normalised"
    ):
        # Paths
        self.main_path = main_path
        self.eigenallele = f"{main_path}.eigenvec.allele"
        self.eigenval_file = f"{main_path}.eigenval"

        self.plink_exe = plink_exe
        self.plink2_exe = plink2_exe

        self.output_dir = output_dir
        self.normalised_dir = normalised_dir
        # Create output directories if they don't exist
        os.makedirs(self.output_dir, exist_ok=True)
        os.makedirs(self.normalised_dir, exist_ok=True)

    # Get the overlapping SNPs
    def overlap_snps(self):
        print("Extracting Reference SNPs -->")
        # Read the eigenvec.allele file and save the second column to a new file
        ref = pd.read_csv(self.eigenallele, delim_whitespace=True)
        # Save the second column to a new file
        ref.iloc[:, 1].to_csv("Reference_SNPs.txt", header=False, index=False)
        print(f"Saved {len(ref)} reference SNPs.\n")

    # Clean and fit PCA
    def clean_fit_pca(self, dtc_file):
        # Pahts
        base_name = os.path.splitext(os.path.basename(dtc_file))[0]
        raw_prefix = f"dtc_raw_{base_name}"
        clean_prefix = f"dtc_clean_{base_name}"

        print(f"Processing: {base_name}")

        # Convert 23andMe to PLINK
        # plink --23file genome.txt --make-bed --out dtc_raw
        cmd_convert = [
            self.plink_exe,
            "--23file", dtc_file,
            "--make-bed",
            "--out", raw_prefix
        ]
        subprocess.run(cmd_convert, check=True)

        # Remove duplicates
        bim = pd.read_csv(
            f"{raw_prefix}.bim",
            sep="\t",
            header=None,
            names=["chr", "snp", "cm", "pos", "a1", "a2"]
        )
        # Identify SNPs where a1 == a2 and also duplicate SNP IDs
        bad_alleles = bim[bim["a1"] == bim["a2"]]["snp"]
        duplicate_ids = bim[bim.duplicated(subset=["snp"])]["snp"]
        # Combine and save bad SNPs
        bad_snps = pd.concat([bad_alleles, duplicate_ids]).drop_duplicates()
        exclude_file = f"exclude_{base_name}.txt"
        bad_snps.to_csv(exclude_file, index=False, header=False)

        # Filter overlap
        # plink --bfile dtc_raw --extract Reference_SNPs.txt --exclude exclude.txt --make-bed --out dtc_clean
        cmd_filter = [
            self.plink_exe,
            "--bfile", raw_prefix,
            "--extract", "Reference_SNPs.txt",
            "--exclude", exclude_file,
            "--make-bed",
            "--out", clean_prefix
        ]
        try:
            subprocess.run(cmd_filter, check=True, stdout=subprocess.DEVNULL)
            print("Filtered overlap and removed bad SNPs.")
        except subprocess.CalledProcessError as e:
            print(f"Error filtering {raw_prefix}: {e}")
            return

        # PCA projection
        n_pc = 15
        start = 7
        end = start + n_pc - 1

        pca_output = os.path.join(self.output_dir, f"{base_name}_pcs")
        # plink2 --bfile dtc_clean --score ref_pca.eigenvec.allele 2 6 header no-mean-imputation --score-col-nums 7-21 --out dtc_pcs/dtc_clean_pcs
        cmd_pca = [
            self.plink2_exe,
            "--bfile", clean_prefix,
            "--score", self.eigenallele, "2", "6",
            "header", "no-mean-imputation",
            "--score-col-nums", f"{start}-{end}",
            "--out", pca_output
        ]

        try:
            subprocess.run(cmd_pca, check=True, stdout=subprocess.DEVNULL)
            print(f"Completed!! Saved to: {pca_output}.sscore\n")
        except subprocess.CalledProcessError as e:
            print(f"Error for {clean_prefix} due to: {e}")

        print(f"PCA complete → {pca_output}.sscore\n")

    # Normalise the PCs
    def normalise_pcs(self):
        # Get the Eigenvalues and the sscore files
        eig = np.loadtxt(self.eigenval_file)
        sscore_files = glob.glob(os.path.join(self.output_dir, "*.sscore"))

        for file in sscore_files:

            df = pd.read_csv(file, sep=r"\s+")

            for i in range(len(eig)):
                score_col = f"SCORE{i+1}_AVG"
                # Normalise by dividing by the square root of the corresponding Eigenvalue
                if score_col in df.columns:
                    df[f"PC{i+1}"] = df[score_col] / np.sqrt(eig[i])
            # Keep only the FID, IID, and the new PC columns
            keep_cols = ["#FID", "IID"] + [f"PC{i+1}" for i in range(len(eig))]
            keep_cols = [c for c in keep_cols if c in df.columns]

            out_path = os.path.join(
                self.normalised_dir,
                os.path.basename(file).replace(".sscore", "_normalised.csv")
            )

            df[keep_cols].to_csv(out_path, index=False)
            print(f"Saved → {out_path}")

        print("All PCA files normalised.")

    # Run the full pipeline
    def run_full_pipeline(self, dtc_file):
        # Get the SNPs list
        self.overlap_snps()
        # Clean & fit PCA
        self.clean_fit_pca(dtc_file)
        # Normalize PCs
        self.normalise_pcs()
        print("Full processing pipeline completed!")

# Class for computing the ancestry and distances
class Ancestry:
    # initial the paths and create output directories
    def __init__(self, user_name=None):
        self.user_name = user_name
        self.df, self.pc_cols, self.target, self.status = self.loading_data()

    def loading_data(self):
        # Load the samples data
        df = pd.read_csv("final_out.csv")
        # Get rid of samples with Ignore in the label
        df = df[~df["Clean_FID"].str.startswith("Ignore", na=False)]
        # Get the PC columns
        pc_cols = [c for c in df.columns if c.startswith("PC")]
        
        # Get the target sample
        if self.user_name:
            pattern = f"normalised/genome_{self.user_name}_pcs_normalised.csv"
        else:
            pattern = "normalised/genome_*_pcs_normalised.csv"
        
        dtc_paths = glob.glob(pattern)
        
        if not dtc_paths:
            # If the user hasn't uploaded a file yet rerurn a message
            return df, pc_cols, None, "Upload file"
            
        dtc_df = pd.read_csv(dtc_paths[0])
        # Get the PC columns
        dtc_pc_cols = [c for c in dtc_df.columns if c.startswith("PC")]
        target = dtc_df[dtc_pc_cols].iloc[0].astype(float).values
        return df, pc_cols, target, "Ready"

    def Compute_nnls(self, df_subset, label):
        # Stop if there is no target file
        if self.target is None:
            return {"label": label, "mixture": "Upload file","raw_components": [], "error": 0.0}
        if len(df_subset) == 0:
            return {"label": label, "mixture": "No samples available", "raw_components": [],"error": 0.0}

        X = df_subset[self.pc_cols].values
        labels = df_subset["Clean_FID"].values

        # Get only the closest samples to save memory
        dist = np.linalg.norm(X - self.target, axis=1)
        N = min(500, len(X))
        idx = np.argsort(dist)[:N]

        X_small = X[idx]
        labels_small = labels[idx]

        # Get the nnls
        weights, _ = scipy_nnls(X_small.T, self.target)

        # incase there is no solution where all weights are 0
        if weights.sum() == 0:
            print(f"\n{label}: No valid solution")
            return {"label": label, "mixture": "No valid solution", "error": 0.0}
            
        # Normalising weights to make sure they add up to 1
        weights /= weights.sum()

        # Choose only the closest 10 samples
        top_k = 10
        top_idx = np.argsort(weights)[-top_k:]
        # Sort the weights in descending order
        results = sorted(
            [(weights[i], labels_small[i]) for i in top_idx],
            reverse=True
        )
        # Filter out samples with small weights. They will show up as 0% when the results are displayed to some decimal place
        threshold = 0.0001  
        results = [
            (w, x) for w, x in results
            if w > threshold
        ]
        # Format into displayable form
        raw_components = [{"weight": w, "label": x} for w, x in results]
        mixture = " + ".join([f"{w*100:.2f}% {x}" for w, x in results])

        # Get the error
        pred = weights @ X_small
        error = np.linalg.norm(pred - self.target)

        return {
            "label": label,
            "mixture": mixture,
            "raw_components": raw_components,
            "error": error
        }

    def compute_distance(self, df_subset, label):
        if self.target is None:
            return {"label": label, "top_closest": []}
        if len(df_subset) == 0:
            return {"label": label, "top_closest": []}

        X = df_subset[self.pc_cols].values
        labels = df_subset["Clean_FID"].values

        # Calculate the Euclidean distance
        dist = np.linalg.norm(X - self.target, axis=1)

        # Get only the top 10 closest
        top_n = min(10, len(dist))
        idx = np.argsort(dist)[:top_n]
        results = [{"label": labels[i], "distance": float(dist[i])} for i in idx]
        return {
            "label": label,
            "top_closest": results
        }

    def ancestry_result(self):
        if self.status != "Ready":
            return None
            
        # Split the data into modern and ancient
        df_modern = self.df[self.df["Broad_Period"] == "Historical/Modern"]
        df_ancient = self.df[self.df["Broad_Period"] != "Historical/Modern"]
        
        # NNLS
        modern_nnls = self.Compute_nnls(df_modern, "Historical/Modern")
        ancient_nnls = self.Compute_nnls(df_ancient, "Ancient")
        
        # Distance
        modern_dist = self.compute_distance(df_modern, "Historical/Modern")
        ancient_dist = self.compute_distance(df_ancient, "Ancient")

        results = {
            "NNLS": {
                "Modern": modern_nnls,
                "Ancient": ancient_nnls
            },
            "Distances": {
                "Modern": modern_dist,
                "Ancient": ancient_dist
            }
        }
        return results

# Class for matching with the predefined traits
class Traits:
    # Load the Json file which contains the SNPs
    def __init__(self, file_path, json_path="Phenotype_SNPs.json"):
        self.file_path = file_path
        self.json_path = json_path
        self.user_df = None
        self.pheno_data = None
    # Get the user file and clean it to get the rsid and genotype columns. The file should be in the format of 23andMe raw data with columns: rsid, chromosome, position, genotype. The function will return a dataframe.
    def load_user_file(self):
	# Get the users uploaded file
        if not os.path.exists(self.file_path):
            raise FileNotFoundError(f"User file not found at {self.file_path}")

        # Clean the file
        self.user_df = pd.read_csv(
            self.file_path,
            sep=r'\s+', 
            comment='#',
            header=None,
            names=['rsid', 'chromosome', 'position', 'genotype'],
            engine='python',
            dtype={'rsid': str, 'genotype': str}
        )
        self.user_df['rsid'] = self.user_df['rsid'].str.strip()
        self.user_df['genotype'] = self.user_df['genotype'].str.strip().str.upper()
        
        return self.user_df
    # Load the JSON file which contains the predefined traits and their associated SNPs. The function will iterate through each trait and check for overlaps with the user's genotype data. It will then add the results to the JSON and return it.
    def load_phenotype_snps(self):
        if self.user_df is None:
            self.load_user_file()
        if not os.path.exists(self.json_path):
            raise FileNotFoundError(f"JSON file not found at {self.json_path}")

        with open(self.json_path, 'r') as f:
            self.pheno_data = json.load(f)
                        # Iterate through each tree in the category

        for x in self.pheno_data:
            category_name = x['name'].strip()
            for y in x.get('children', []):
                # Extract the trees from the JSOON
                try:
                    rsid = y['name'].strip()
                    gene = y['children'][0]
                    trait = gene['children'][0]
                    allele = trait['children'][0]

                    allele_name = allele['name'].strip().upper()
                    gene_name = gene['name'].strip()
                    trait_name = trait['name'].strip()
                except AttributeError:
                    continue
                # Check overlap
                match = self.user_df[self.user_df['rsid'] == rsid]
                # Matching logic
                if match.empty:
                    status = -1
                    user_gt = "--"
                else:
                    # count the number of matches
                    user_gt = match.iloc[0]['genotype']
                    count = user_gt.count(allele_name)
                    if user_gt in ["--", "00", "II", "DD"]:
                        status = -1
                    else:
                        status = count
                y['results']={
                    'Region': category_name,
                    'RsID': rsid,
                    'gene': gene_name,
                    'trait': trait_name,
                    'target': allele_name,
                    'user_gt': user_gt,
                    'status': status
                }
        return self.pheno_data  
    # Do the same for the health traits
    def get_health_results(self):
        try:
            # Load the  health JSON
            with open('health_markers.json', 'r') as f:
                health_data = json.load(f)

            # Iterate through each category
            for category in health_data:
                category_name = category.get('name', 'Unknown').strip()
                
                # Iterate through each tree in the category
                for x in category.get('children', []):
                    try:
                        rsid = x['name'].strip()
                        
                        gene_node = x['children'][0]
                        trait_node = gene_node['children'][0]
                        allele_node = trait_node['children'][0]

                        # Extract names for display and matching
                        gene_name = gene_node['name'].strip()
                        trait_name = trait_node['name'].strip()
                        allele_name = allele_node['name'].strip().upper()
                    except (KeyError, IndexError, AttributeError):
                        continue

                    # Check overlap
                    match = self.user_df[self.user_df['rsid'] == rsid]
                    # Matching logic
                    if match.empty:
                        status = -1
                        user_gt = "--"
                    else:
                        user_gt = match.iloc[0]['genotype'].upper()
                        
                        if user_gt in ["--", "00", "II", "DD", "??"]:
                            status = -1
                        else:
                            # Count matches: 0 No Match, 1 partial and 2 full
                            status = user_gt.count(allele_name)

                    
                    x['results'] = {
                        'Region': category_name,
                        'RsID': rsid,
                        'gene': gene_name,
                        'trait': trait_name,
                        'target': allele_name,
                        'user_gt': user_gt,
                        'status': status
                    }

            return health_data

        except FileNotFoundError:
            print("Error: health_markers.json not found.")
            return []
        except Exception as e:
            print(f"Error processing results: {e}")
            return []

