import streamlit as st
from ui_utils import (
    render_instructions,
    render_load_data,
    render_visualization,
    render_all_datasets,
    render_ai_assistant
)

st.set_page_config(page_title="RNA-seq Viewer", layout="wide")

st.sidebar.title("Navigation")
menu = st.sidebar.radio("Select a section:", [
    "📘 Instructions", "📁 Load Data", "📊 Visualize", "🗂️ All Datasets", "🤖 AI Assistant"
])

if menu == "📘 Instructions":
    render_instructions()

elif menu == "📁 Load Data":
    render_load_data()

elif menu == "📊 Visualize":
    render_visualization()

elif menu == "🤖 AI Assistant":
    render_ai_assistant()

st.markdown("---")
st.info("Let's shine in this hackathon, team Goat! 🐐")
