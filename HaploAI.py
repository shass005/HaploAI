import os
import streamlit as st
from streamlit_option_menu import option_menu
import pandas as pd
#from st_aggrid import AgGrid, GridOptionsBuilder, JsCode
#from st_aggrid.shared import GridUpdateMode
import streamlit.components.v1 as components
import numpy as np
from Target import TargetProcessor
from Target import Ancestry

tarp = TargetProcessor()
st.set_page_config(page_title="HaploAI Dashboard", layout="wide")

with st.sidebar:
    light_mode = st.toggle("Light Mode")
    #if light_mode :
     #   st.write("Coming Soon!")
    
    selected = option_menu(
        menu_title = "Main Menu",
        options = ["Home","Upload","Ancestry", "Phenotype", "Health"],
        icons = ["house","bi-cloud-arrow-down", "person-vcard", "eye","lungs"],
        menu_icon = "body-text",
        default_index = 0,
    )

    if light_mode:
        st.markdown("""
            <style>
            /* Light mode styling */
            body, .main, .block-container {
                background-color: #f0f0f0 !important;
                color: black !important;
            }

            .stApp {
                background-color: #f0f0f0 !important;
                color: black !important;
            }
            h1{
                color: black !important;
            }
            .st-emotion-cache-15hul6a {
                border-color:black !important;
            }
            /*sidebar*/
            .st-emotion-cache-155jwzh{
                background-color: #1e1e1e !important;
            }
            .st-emotion-cache-1gkg1x5 p, .st-emotion-cache-1gkg1x5 ol, .st-emotion-cache-1gkg1x5 ul, .st-emotion-cache-1gkg1x5 dl, .st-emotion-cache-1gkg1x5 li{
                color: white !important;
            }
            p{
                color: black !important;
            }  
            .st-emotion-cache-fm8pe0 p {
                color: white !important;
            }
            .st-emotion-cache-13na8ym {
                border-color:black !important;        
            }
            .card {
                background-color:white !important;
            }
            .st-emotion-cache-15hul6a {
                background-color:white !important;
            }
            .st-emotion-cache-fm8pe0 p {
                color:#143ab2 !important;
            }
            /* Sidebar light */
            section[data-testid="stSidebar"] {
                background-color: #f8f9fa !important;
                color: black !important;
            }

            .nav-link, .nav-link i {
                color: black !important;
            }
            .nav-link.active {
                background-color: #e0e0e0 !important;
                font-weight: bold;
            }

            </style>
        """, unsafe_allow_html=True)
    else:
        st.markdown("""
            <style>
            /* Dark mode styling */
            body, .main, .block-container {
                background-color: #1e1e1e !important;
                color: white !important;
            }

            .stApp {
                background-color: #1e1e1e !important;
                color: white !important;
            }

            /* Sidebar dark */
            section[data-testid="stSidebar"] {
                background-color: #1e1e1e !important;
                color: white !important;
            }

            .nav-link, .nav-link i {
                color: white !important;
            }
            .nav-link.active {
                background-color: #333333 !important;
                font-weight: bold;
            }
            </style>
        """, unsafe_allow_html=True)

if selected == "Home":
    st.title("Home")
    st.write("Welcome to HaploAI, .... To use the service upload your file in the upload tab and check the repective tabs for the: Ancestry, Phenotypical and Health results. ")    

if selected == "Upload":
    st.title("Upload")
    # Get the users name and store it
    user_name = st.text_input("Enter your name:", value=st.session_state.get('user_input_name', ""))
    uploaded_file = st.file_uploader("Upload your 23andMe raw DNA file", type=["txt"])

    # Ensure both file and name have inputs
    if uploaded_file is not None and user_name and user_name.strip():
            # Store the file
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
                    st.success(f"Analysis complete for {user_name}!")
            else:
                st.info(f"Data for {user_name} is already processed and ready.")

if selected == "Ancestry":
    user = st.session_state.get('processed_user')
    
    if not user:
        st.warning("DNA Data Required\n Please upload and process your file in the **Upload** tab first.")
    else:
        # Custom css for the ancestry resutls
        st.markdown("""
            <style>
            /* Modern Section  */
            .modern-card {
                background-color: rgba(0, 104, 201, 0.1);
                border-left: 5px solid #0068c9;
                padding: 15px;
                border-radius: 10px;
                margin-bottom: 10px;
            }
            /* Ancient Section  */
            .ancient-card {
                background-color: rgba(255, 166, 0, 0.1);
                border-left: 5px solid #ffa600;
                padding: 15px;
                border-radius: 10px;
                margin-bottom: 10px;
            }
            /* Distance Cards */
            .dist-card {
                background: #262730;
                border: 1px solid #464855;
                border-radius: 8px;
                padding: 10px;
                text-align: center;
                margin: 5px;
            }
            </style>
        """, unsafe_allow_html=True)

        with st.spinner("Please Wait, generating Genetic Dashboard..."):
            analyzer = Ancestry(user_name=user)
            results = analyzer.ancestry_result()

        if results:
            st.title(f"Genomic Dashboard: {user}")
            
            # This is top half, ancestry
            col_top_left, col_top_right = st.columns(2, gap="medium")

            # Moderls results
            with col_top_left:
                st.markdown("🏙️ Modern Populations")
                # Display error results
                st.metric("Model Fit Error", f"{results['NNLS']['Modern']['error']:.5f}")
                with st.container(border=True):
                    for item in results['NNLS']['Modern']['raw_components']:
                        pct = item['weight'] * 100
                        st.markdown(f"**{item['label']}** ({pct:.1f}%)")
                        st.progress(item['weight'])

            # Ancient results
            with col_top_right:
                st.markdown("🏛️ Ancient Populations")
                st.metric("Model Fit Error", f"{results['NNLS']['Ancient']['error']:.5f}")
                with st.container(border=True):
                    for item in results['NNLS']['Ancient']['raw_components']:
                        pct = item['weight'] * 100
                        st.markdown(f"**{item['label']}** ({pct:.1f}%)")
                        st.progress(item['weight'])

            # Seperate the page for clean dashboard look
            st.markdown("---")
            st.header("📍 Genetic Distance")
            st.write("Lower distances indicate genetic proximity to the given population.")

            # Lower half, distances
            col_bot_left, col_bot_right = st.columns(2, gap="medium")

            with col_bot_left:
                st.subheader("Closest Modern Matches")
               # Display the results in mini cards with streamlit html injections
                grid_m = st.columns(2)
                for i, r in enumerate(results["Distances"]["Modern"]["top_closest"]):
                    with grid_m[i % 2]: 
                        st.markdown(f"""
                            <div class="modern-card">
                                <small>Rank #{i+1}</small><br>
                                <strong>{r['label']}</strong><br>
                                <code style="color:#0068c9;">{r['distance']:.5f}</code>
                            </div>
                        """, unsafe_allow_html=True)

            with col_bot_right:
                st.subheader("Closest Ancient Matches")
                
                grid_a = st.columns(2)
                for i, r in enumerate(results["Distances"]["Ancient"]["top_closest"]):
                    with grid_a[i % 2]:
                        st.markdown(f"""
                            <div class="ancient-card">
                                <small>Rank #{i+1}</small><br>
                                <strong>{r['label']}</strong><br>
                                <code style="color:#ffa600;">{r['distance']:.5f}</code>
                            </div>
                        """, unsafe_allow_html=True)
                        
        else:
            st.error("Data could not be loaded. Please ensure processing completed successfully.")       

if selected == "Phenotype":
    user = st.session_state.get('processed_user')
    st.title(f"Phenotype Results: {user}")
