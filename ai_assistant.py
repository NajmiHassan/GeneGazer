import streamlit as st
import google.generativeai as genai
import os

class GeminiAssistant:
    """Gemini AI Assistant for RNA-seq analysis"""
    
    def __init__(self):
        self.model = None
        self.setup_gemini()
    
    def setup_gemini(self):
        """Setup Gemini AI with API key"""
        api_key = st.secrets.get("GEMINI_API_KEY") or os.getenv("GEMINI_API_KEY")
        if not api_key:
            st.error("⚠️ Gemini API key not found. Please add GEMINI_API_KEY to your Streamlit secrets or environment variables.")
            return None
        
        try:
            genai.configure(api_key=api_key)
            self.model = genai.GenerativeModel('gemini-1.5-flash')
            return self.model
        except Exception as e:
            st.error(f"Error setting up Gemini: {str(e)}")
            return None
    
    def get_dataset_context(self):
        """Get context about currently loaded dataset"""
        if 'adata' not in st.session_state:
            return "No dataset is currently loaded."
        
        adata = st.session_state['adata']
        context = f"""
        Current dataset information:
        - Number of cells: {adata.n_obs}
        - Number of genes: {adata.n_vars}
        - Available clusters: {len(adata.obs['leiden'].unique()) if 'leiden' in adata.obs else 'Not clustered yet'}
        - Top variable genes: {', '.join(adata.var_names[:10].tolist())}
        """
        return context
    
    def generate_response(self, user_input):
        """Generate AI response with context"""
        if self.model is None:
            return "AI Assistant is not available. Please check your API key configuration."
        
        try:
            # Create context-aware prompt
            dataset_context = self.get_dataset_context()
            
            system_prompt = f"""You are an expert bioinformatics assistant specializing in single-cell RNA sequencing analysis. 
            You help researchers understand their data, interpret results, and provide biological insights.
            
            Current user's dataset context:
            {dataset_context}
            
            Please provide helpful, accurate, and scientifically sound responses. If discussing specific genes or pathways, 
            provide biological context. If discussing analysis methods, explain the rationale behind different approaches.
            
            User question: {user_input}"""
            
            # Generate response
            response = self.model.generate_content(system_prompt)
            return response.text
            
        except Exception as e:
            return f"Sorry, I encountered an error: {str(e)}"
    
    def initialize_chat_history(self):
        """Initialize chat history if it doesn't exist"""
        if "chat_history" not in st.session_state:
            st.session_state.chat_history = []
    
    def display_greeting(self):
        """Display initial greeting message"""
        greeting = """👋 Hi! I'm your AI assistant powered by Gemini. I can help you with:

📊 **Data Analysis**: Interpret your single-cell RNA-seq results
🧬 **Biology**: Explain gene functions, pathways, and cell types  
🔬 **Methods**: Discuss preprocessing, clustering, and visualization techniques
📈 **Statistics**: Help understand your data metrics and quality control

Feel free to ask me anything about your RNA-seq data or bioinformatics in general!"""
        
        with st.chat_message("assistant"):
            st.markdown(greeting)
        
        st.session_state.chat_history.append({"role": "assistant", "content": greeting})
    
    def display_chat_history(self):
        """Display existing chat history"""
        for message in st.session_state.chat_history:
            with st.chat_message(message["role"]):
                st.markdown(message["content"])
    
    def add_message_to_history(self, role, content):
        """Add message to chat history"""
        st.session_state.chat_history.append({"role": role, "content": content})
    
    def clear_chat_history(self):
        """Clear chat history"""
        st.session_state.chat_history = []
    
    def render_ai_assistant_tab(self):
        """Render the complete AI Assistant tab"""
        st.title("🧬 Gemini AI Assistant for RNA-seq")
        
        # Check if model is available
        if self.model is None:
            st.stop()
        
        # Initialize chat history
        self.initialize_chat_history()
        
        # Display chat history
        self.display_chat_history()
        
        # Initial greeting if no chat history
        if not st.session_state.chat_history:
            self.display_greeting()
        
        # Chat input
        if user_input := st.chat_input("Ask me anything about RNA-seq..."):
            # Add user message to chat history
            self.add_message_to_history("user", user_input)
            
            # Display user message
            with st.chat_message("user"):
                st.markdown(user_input)
            
            # Generate AI response
            with st.chat_message("assistant"):
                with st.spinner("🤔 Thinking..."):
                    ai_response = self.generate_response(user_input)
                    
                    # Display response
                    if ai_response.startswith("Sorry, I encountered an error"):
                        st.error(ai_response)
                    else:
                        st.markdown(ai_response)
                    
                    # Add to chat history
                    self.add_message_to_history("assistant", ai_response)
        
        # Clear chat button
        if st.button("🗑️ Clear Chat History"):
            self.clear_chat_history()
            st.rerun()