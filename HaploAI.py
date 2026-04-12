import os
import streamlit as st
from streamlit_option_menu import option_menu
import pandas as pd
import streamlit.components.v1 as components
import numpy as np
import plotly.express as px 
from Target import TargetProcessor
from Target import Ancestry
from Target import Traits

tarp = TargetProcessor()
st.set_page_config(page_title="HaploAI Dashboard", layout="wide", initial_sidebar_state="expanded")

with st.sidebar:
    light_mode = st.toggle("Light Mode")
    # Main menu with icons    
    selected = option_menu(
        menu_title = "Main Menu",
        options = ["Home", "Upload", "Dashboard"],
        icons = ["house", "bi-cloud-arrow-down", "columns-gap"],
        menu_icon = "body-text",
        default_index = 0,
    )

    if light_mode:
        st.markdown("""
            <style>
            /* Light mode styling */
            body, .main, .block-container { background-color: #f0f0f0 !important; color: black !important; }
            .stApp { background-color: #f0f0f0 !important; color: black !important; }
            h1, p { color: black !important; }
            .st-emotion-cache-15hul6a { border-color:black !important; background-color:white !important; }
            .st-emotion-cache-155jwzh { background-color: #1e1e1e !important; }
            .st-emotion-cache-1gkg1x5 p, .st-emotion-cache-1gkg1x5 ol, .st-emotion-cache-1gkg1x5 ul, .st-emotion-cache-1gkg1x5 dl, .st-emotion-cache-1gkg1x5 li { color: white !important; }
            .st-emotion-cache-fm8pe0 p { color:#143ab2 !important; }
            section[data-testid="stSidebar"] { background-color: #f8f9fa !important; color: black !important; }
            .nav-link, .nav-link i { color: black !important; }
            .nav-link.active { background-color: #e0e0e0 !important; font-weight: bold; }
            </style>
        """, unsafe_allow_html=True)
    else:
        st.markdown("""
            <style>
            /* Dark mode styling */
            body, .main, .block-container { background-color: #1e1e1e !important; color: white !important; }
            .stApp { background-color: #1e1e1e !important; color: white !important; }
            section[data-testid="stSidebar"] { background-color: #1e1e1e !important; color: white !important; }
            .nav-link, .nav-link i { color: white !important; }
            .nav-link.active { background-color: #333333 !important; font-weight: bold; }
            
            /* Category Header Layout */
            .cat-header { display: flex; justify-content: space-between; align-items: center; background: #262730; padding: 15px; border-radius: 10px 10px 0 0; border-bottom: 2px solid #333; }
            .cat-title { font-size: 1.5rem; font-weight: bold; color: white; margin: 0; }
            .cat-prediction { font-size: 0.9rem; color: #aaa; margin: 0; }
            .color-bar { width: 80px; height: 12px; border-radius: 6px; }
            
            /* SNP Cards */
            .snp-card { padding: 12px; margin: 8px 0; border-radius: 8px; border-left: 5px solid; }
            .status-2 { background: rgba(40, 167, 69, 0.15); border-left-color: #28a745; }
            .status-1 { background: rgba(255, 193, 7, 0.15); border-left-color: #ffc107; }
            .status-0 { background: rgba(220, 53, 69, 0.15); border-left-color: #dc3545; }
            .status--1 { background: rgba(108, 117, 125, 0.15); border-left-color: #6c757d; }
            
            /* Distance Cards */
            .modern-card { background-color: rgba(0, 104, 201, 0.1); border-left: 5px solid #0068c9; padding: 15px; border-radius: 10px; margin-bottom: 10px; }
            .ancient-card { background-color: rgba(255, 166, 0, 0.1); border-left: 5px solid #ffa600; padding: 15px; border-radius: 10px; margin-bottom: 10px; }
            .profile-card { background: #262730; border: 1px solid #464855; border-radius: 10px; padding: 20px; height: 100%; }
            </style>
        """, unsafe_allow_html=True)
# Home page
if selected == "Home":
    st.title("Home")
    st.write("Welcome to HaploAI..... To use the service, upload your file in the Upload tab and check the Dashboard for your results.")    
# Upload page
elif selected == "Upload":
    st.title("Upload")
    # Get the users name and store it
    user_name = st.text_input("Enter your name:", value=st.session_state.get('user_input_name', ""))
    uploaded_file = st.file_uploader("Upload your 23andMe raw DNA file", type=["txt"])
    # Ensure user has entered both a file and their name
    if uploaded_file is not None and user_name and user_name.strip():
        st.session_state['user_input_name'] = user_name.strip()
        name = f"genome_{user_name.strip()}"
        os.makedirs("23andme", exist_ok=True)
        dtc_file_path = f"23andme/{name}.txt"
        
        if st.session_state.get('processed_user') != user_name.strip():
            with open(dtc_file_path, "wb") as f:
                f.write(uploaded_file.getbuffer())

            st.success("File uploaded successfully!")
            with st.spinner("Processing DNA data..."):
                tarp.run_full_pipeline(dtc_file_path)
                st.session_state['processed_user'] = user_name.strip()
                st.session_state['user_file_path'] = dtc_file_path
                st.success(f"Analysis complete for {user_name}! Navigate to the Dashboard.")
        else:
            st.info(f"Data for {user_name} is already processed and ready.")
# Dashboard page for results
elif selected == "Dashboard":
    # Get the user and their file
    user = st.session_state.get('processed_user')
    user_file = st.session_state.get('user_file_path')
    # In case the user navigates to the results page first
    # Display warning prompt
    if not user or not user_file:
        st.warning("DNA Data Required\nPlease upload and process your file in the **Upload** tab first.")
    else:
        st.title(f"Genomic Dashboard: {user}")
        
        with st.spinner("Loading Dashboard Data..."):
            # Load Data
            analyzer = Ancestry(user_name=user)
            anc_results = analyzer.ancestry_result()
            
            engine = Traits(file_path=user_file)
            pheno_results = engine.load_phenotype_snps()

        if anc_results:
           # First row
           # Display the profile and ancestry results
            top_left, top_mid, top_right = st.columns([2, 2, 1], gap="large")
            # Modenr ancestry section
            with top_left:
                st.markdown("⛫ Modern Ancestry")
                modern_data = anc_results['NNLS']['Modern']['raw_components']
                # If the calculation produces some meaningful results for the user
                # Plot the results as a pie chart
                if modern_data:
                    df_mod = pd.DataFrame(modern_data)
                    fig_mod = px.pie(df_mod, values='weight', names='label', hole=0.4)
                    fig_mod.update_layout(margin=dict(t=20, b=20, l=0, r=0), showlegend=True)
                    st.plotly_chart(fig_mod, use_container_width=True)
                else:
                    st.write("Not enough data!")
            # Ancient ancestry results
            with top_mid:
                st.markdown("𓉱 Ancient Ancestry")
                ancient_data = anc_results['NNLS']['Ancient']['raw_components']
                if ancient_data:
                    df_anc = pd.DataFrame(ancient_data)
                    fig_anc = px.pie(df_anc, values='weight', names='label', hole=0.4)
                    fig_anc.update_layout(margin=dict(t=20, b=20, l=0, r=0), showlegend=True)
                    st.plotly_chart(fig_anc, use_container_width=True)
                else:
                    st.write("Not enough data!")

            with top_right:
                st.markdown("Profile")
                st.markdown(f"""
                    <div class="profile-card">
                        <p style="font-size: 1.1rem;"><strong>Name:</strong> <br><span style="color:#0068c9;">{user}</span></p>
                        <p style="font-size: 1.1rem;"><strong>YDNA:</strong> <br><span style="color:#0068c9;">xxx</span></p>
                        <p style="font-size: 1.1rem;"><strong>mtDNA:</strong> <br><span style="color:#0068c9;">xxx</span></p>
                    </div>
                """, unsafe_allow_html=True)
            
            st.markdown("---")

           # Second row
           # Distances
            st.markdown("⚲ Genetic Distances")
            dist_left, dist_right = st.columns(2, gap="large")
            # Modern section
            with dist_left:
                st.markdown("Modern Distances")
                with st.container(border=True):
                    grid_m = st.columns(2)
                    # The page doesn't look good with too many distance results so reduce to top 5 only
                    for i, r in enumerate(anc_results["Distances"]["Modern"]["top_closest"][:5]): 
                        with grid_m[i % 2]: 
                            st.markdown(f"""
                                <div class="modern-card">
                                    <small>Rank #{i+1}</small><br>
                                    <strong>{r['label']}</strong><br>
                                    <code style="color:#0068c9;">{r['distance']:.5f}</code>
                                </div>
                            """, unsafe_allow_html=True)

            with dist_right:
                st.markdown("Ancient Distances")
                with st.container(border=True):
                    grid_a = st.columns(2)
                    for i, r in enumerate(anc_results["Distances"]["Ancient"]["top_closest"][:5]):
                        with grid_a[i % 2]:
                            st.markdown(f"""
                                <div class="ancient-card">
                                    <small>Rank #{i+1}</small><br>
                                    <strong>{r['label']}</strong><br>
                                    <code style="color:#ffa600;">{r['distance']:.5f}</code>
                                </div>
                            """, unsafe_allow_html=True)
            
            st.markdown("---")

            # Third part for phenotypical results
            st.markdown("👁 Phenotypical Traits")
            
            # Create divs to house the results for the three regions
            display_categories = ["Hair", "Skin", "Eye"]
            pheno_cols = st.columns(3, gap="large")

            for i, cat_name in enumerate(display_categories):
                category_data = next((c for c in pheno_results if c['name'] == cat_name), None)
                
                with pheno_cols[i]:
                    with st.container(border=True):
                        if category_data and category_data.get('children'):
                            best_match = max(category_data['children'], key=lambda x: x['results']['status'])
                            prediction = best_match['results']['trait'] if best_match['results']['status'] > 0 else "Indeterminate"
                            bar_color = "#28a745" if best_match['results']['status'] == 2 else "#ffc107" if best_match['results']['status'] == 1 else "#6c757d"

                            st.markdown(f"""
                                <div class="cat-header">
                                    <div>
                                        <p class="cat-title">{cat_name}</p>
                                        <p class="cat-prediction">Likely: {prediction}</p>
                                    </div>
                                    <div class="color-bar" style="background-color: {bar_color};"></div>
                                </div>
                            """, unsafe_allow_html=True)

                            for snp in category_data['children']:
                                res = snp['results']
                                status_cls = f"status-{res['status']}"
                                st_text = "Strong Match" if res['status'] == 2 else "Partial Match" if res['status'] == 1 else "No Match" if res['status'] == 0 else "SNP Absent"

                                st.markdown(f"""
                                    <div class="snp-card {status_cls}">
                                        <div style="display:flex; justify-content:space-between; font-size:0.75rem; opacity:0.8;">
                                            <span>{res['RsID']}</span>
                                            <span>{res['gene']}</span>
                                        </div>
                                        <div style="font-weight:bold; font-size:0.9rem; margin: 4px 0;">{res['trait']}</div>
                                        <div style="font-size:0.7rem;">Target: {res['target']} | User: {res['user_gt']}</div>
                                        <div style="font-size:0.7rem; font-style:italic; margin-top:5px;">{st_text}</div>
                                    </div>
                                """, unsafe_allow_html=True)
                        else:
                            st.markdown(f"""
                                <div class="cat-header">
                                    <div><p class="cat-title">{cat_name}</p></div>
                                </div>
                                <div style="padding: 15px;">No markers found.</div>
                            """, unsafe_allow_html=True)