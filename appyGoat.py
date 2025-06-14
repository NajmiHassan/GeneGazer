import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt
import tempfile

st.set_page_config(page_title="Visualizador RNA-seq", layout="wide")

def carregar_e_preprocessar(adata):
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata = adata[:, adata.var.highly_variable]
    sc.pp.scale(adata, max_value=10)
    return adata

def aplicar_pca_umap_cluster(adata):
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    sc.tl.leiden(adata)
    return adata

st.sidebar.title("Navegação")
menu = st.sidebar.radio("Escolha uma seção:", ["📘 Instruções", "📁 Carregar Dados", "📊 Visualizar", "🤖 Agente de IA"])

if menu == "📘 Instruções":
    st.title("Visualizador de RNA-seq de Célula Única")
    st.markdown("""
    Esta ferramenta permite explorar dados de RNA-seq de célula única de forma simples:

    1. Vá até **Carregar Dados** e envie um arquivo `.h5ad` ou use um dataset de exemplo.
    2. Vá até **Visualizar** para gerar gráficos interativos:
       - Agrupamento (UMAP)
       - Mapa de calor por gene
    3. Use o **Agente de IA** para tirar dúvidas.

    **Formato aceito**:
    Arquivos `.h5ad` exportados pelo Scanpy ou Seurat.
    """)
    st.image("https://i.imgur.com/zVfGZkP.png", caption="Exemplo de visualização UMAP")

elif menu == "📁 Carregar Dados":
    st.title("Carregar Dados de RNA-seq")

    arquivo_enviado = st.file_uploader("Envie um arquivo `.h5ad`", type=["h5ad"])
    if arquivo_enviado:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".h5ad") as tmp:
            tmp.write(arquivo_enviado.read())
            caminho_temp = tmp.name
        try:
            adata = sc.read_h5ad(caminho_temp)
            adata = carregar_e_preprocessar(adata)
            adata = aplicar_pca_umap_cluster(adata)
            st.session_state['adata'] = adata
            st.success("Arquivo carregado e processado com sucesso!")
        except Exception as e:
            st.error(f"Erro durante o processamento: {str(e)}")

    st.markdown("### Ou escolha um dos conjuntos de dados públicos:")
    col1, col2 = st.columns(2)

    with col1:
        if st.button("Usar PBMC 3k (Scanpy)"):
            try:
                adata = sc.datasets.pbmc3k()
                adata = carregar_e_preprocessar(adata)
                adata = aplicar_pca_umap_cluster(adata)
                st.session_state['adata'] = adata
                st.success("PBMC 3k carregado e processado!")
            except Exception as e:
                st.error(f"Erro ao carregar PBMC 3k: {str(e)}")

    with col2:
        if st.button("Usar Cérebro de Camundongo (Scanpy)"):
            try:
                adata = sc.datasets.mouse_brain_cortex_1k()
                adata = carregar_e_preprocessar(adata)
                adata = aplicar_pca_umap_cluster(adata)
                st.session_state['adata'] = adata
                st.success("Cérebro de camundongo carregado e processado!")
            except Exception as e:
                st.error(f"Erro ao carregar cérebro de camundongo: {str(e)}")

elif menu == "📊 Visualizar":
    st.title("Visualizações")

    if 'adata' not in st.session_state:
        st.warning("Por favor, carregue um arquivo na aba 'Carregar Dados'.")
    else:
        adata = st.session_state['adata']
        col1, col2 = st.columns(2)

        with col1:
            st.subheader("Visualização UMAP")
            sc.pl.umap(adata, color='leiden', show=False)
            st.pyplot(plt.gcf())

        with col2:
            st.subheader("Mapa de Calor por Gene")
            gene = st.text_input("Digite o nome de um gene:", "IL7R")
            if gene in adata.var_names:
                sc.pl.matrixplot(adata, var_names=[gene], groupby="leiden", cmap="viridis", use_raw=False, show=False)
                st.pyplot(plt.gcf())
            else:
                st.error("Gene não encontrado.")

elif menu == "🤖 Agente de IA":
    st.title("Agente de IA: Pergunte sobre RNA-seq")

    with st.chat_message("assistant"):
        st.markdown("Olá! Pergunte qualquer coisa sobre RNA-seq, genes ou como usar esta ferramenta!")

    pergunta = st.chat_input("Digite sua pergunta aqui...")
    if pergunta:
        with st.chat_message("user"):
            st.markdown(pergunta)

        with st.chat_message("assistant"):
            if "gene" in pergunta.lower():
                st.markdown("Um gene é uma sequência de DNA que contém instruções para produzir proteínas.")
            elif "umap" in pergunta.lower():
                st.markdown("UMAP é uma técnica que reduz dados complexos para 2D e facilita a visualização.")
            elif "h5ad" in pergunta.lower():
                st.markdown("O formato `.h5ad` é usado pelo Scanpy para armazenar dados de expressão genética.")
            else:
                st.markdown("Sou um agente simples, mas posso ajudar com perguntas básicas sobre RNA-seq.")

st.markdown("---")
st.info("Vamos brilhar nesse hackathon, equipe Goat!")
