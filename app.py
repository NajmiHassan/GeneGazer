import streamlit as st
from ui_utils import render_instructions, render_load_data, render_ai_assistant
from visualization_ui import render_visualization

# Set page config
st.set_page_config(page_title="GeneGazer", layout="wide")

# Project title and tagline
st.title("🧬 GeneGazer: Interactive Single-Cell RNA-seq Explorer")
st.caption("Visualize, explore, and analyze single-cell data effortlessly.")

# Main tab navigation
tab1, tab2, tab3, tab4 = st.tabs([
    "📘 Instructions", 
    "📁 Load Data", 
    "📊 Visualize", 
    "🤖 AI Assistant"
])

with tab1:
    render_instructions()

with tab2:
    render_load_data()

with tab3:
    render_visualization()

with tab4:
    render_ai_assistant()

# Optional footer
st.markdown("---")
st.info("Made by Team GOAT 🧬")