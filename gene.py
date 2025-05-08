import streamlit as st
from Bio.Blast import NCBIWWW
from Bio import SeqIO
import requests
import pubchempy as pcp
from Bio.KEGG import REST

# Set page configuration
st.set_page_config(
    page_title="SynoBio Toolkit",
    page_icon="🧬",
    layout="wide"
)

# Custom CSS for better UI
st.markdown(
    """
    <style>
        :root {
            --primary: #2c3e50;
            --secondary: #3498db;
            --accent: #e74c3c;
            --light: #ecf0f1;
            --dark: #2c3e50;
            --success: #2ecc71;
        }
        .main {
            background: linear-gradient(135deg, #f5f7fa 0%, #e4e8f0 100%);
        }
        .sidebar .sidebar-content {
            background-color: #D8BFD8;
        }
        [data-testid="stHeader"] {
            background-color: var(--primary);
        }
        [data-testid="stToolbar"] {
            display: none;
        }
        .stTabs [data-baseweb="tab-list"] {
            gap: 0;
            padding: 0 1rem;
        }
        .stTabs [data-baseweb="tab"] {
            padding: 12px 24px;
            margin: 0;
            border-radius: 4px 4px 0 0;
            background-color: var(--light);
            transition: all 0.2s;
        }
        .stTabs [aria-selected="true"] {
            background-color: var(--secondary);
            color: white !important;
            font-weight: 600;
        }
        .stTabs [aria-selected="false"] {
            color: var(--dark);
        }
        .sequence-box {
            font-family: monospace;
            background-color: #ffffff;
            padding: 15px;
            border-radius: 5px;
            border: 1px solid #dee2e6;
            margin: 10px 0;
            white-space: pre-wrap;
            word-break: break-all;
        }
        .feature-card {
            background: white;
            border-radius: 8px;
            padding: 20px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            margin-bottom: 20px;
            height: 100%;
            transition: transform 0.3s ease;
        }
        .feature-card:hover {
            transform: translateY(-5px);
            box-shadow: 0 5px 15px rgba(0,0,0,0.1);
        }
        .metric-card {
            background: var(--light);
            border-radius: 8px;
            padding: 15px;
            text-align: center;
            border-left: 4px solid var(--secondary);
        }
        .footer {
            text-align: center;
            padding: 20px;
            color: #7f8c8d;
            font-size: 0.9rem;
            margin-top: 40px;
            border-top: 1px solid #eee;
        }
        .hero {
            background: linear-gradient(135deg, var(--primary) 0%, var(--secondary) 100%);
            color: white;
            padding: 3rem 2rem;
            border-radius: 10px;
            margin-bottom: 2rem;
        }
        .team-card {
            background: white;
            border-radius: 8px;
            padding: 20px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            margin-bottom: 20px;
            text-align: center;
        }
        .team-card img {
            width: 120px;
            height: 120px;
            border-radius: 50%;
            object-fit: cover;
            margin-bottom: 15px;
            border: 3px solid var(--secondary);
        }
        .tooltip {
            position: relative;
            display: inline-block;
            border-bottom: 1px dotted black;
        }
        .tooltip .tooltiptext {
            visibility: hidden;
            width: 200px;
            background-color: #555;
            color: #fff;
            text-align: center;
            border-radius: 6px;
            padding: 5px;
            position: absolute;
            z-index: 1;
            bottom: 125%;
            left: 50%;
            margin-left: -100px;
            opacity: 0;
            transition: opacity 0.3s;
        }
        .tooltip:hover .tooltiptext {
            visibility: visible;
            opacity: 1;
        }
        .badge {
            display: inline-block;
            padding: 0.25em 0.4em;
            font-size: 75%;
            font-weight: 700;
            line-height: 1;
            text-align: center;
            white-space: nowrap;
            vertical-align: baseline;
            border-radius: 0.25rem;
            background-color: var(--accent);
            color: white;
        }
        body {
            font-family: 'Arial', sans-serif;
            color: #333333;
        }
        .title {
            color: #6A5ACD;
            font-size: 40px;
            text-align: center;
            margin-bottom: 20px;
        }
        .tool-header {
            color: #6A5ACD;
            border-bottom: 2px solid #6A5ACD;
            padding-bottom: 10px;
        }
        .info-box {
            background-color: #ffffff;
            padding: 15px;
            border-radius: 10px;
            margin-bottom: 20px;
            box-shadow: 0 2px 5px rgba(0,0,0,0.1);
        }
        .stTextInput>div>div>input {
            background-color: #ffffff;
        }
        .stSelectbox>div>div>select {
            background-color: #ffffff;
        }
        .author-section {
            background-color: #ffffff;
            padding: 20px;
            border-radius: 10px;
            margin: 20px 0;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }
        .server-info {
            background-color: #ffffff;
            padding: 20px;
            border-radius: 10px;
            margin: 20px 0;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }
        .mentorship-section {
            background-color: #ffffff;
            padding: 20px;
            border-radius: 10px;
            margin: 20px 0;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }
    </style>
    """,
    unsafe_allow_html=True,
)

# Cache API requests
@st.cache_data(ttl=3600)
def run_blast(sequence, program):
    try:
        result_handle = NCBIWWW.qblast(program=program, database="nr", sequence=sequence)
        return result_handle.read()
    except Exception as e:
        st.error(f"BLAST failed: {e}")
        return None

@st.cache_data(ttl=3600)
def fetch_uniprot(protein_id):
    try:
        response = requests.get(f"https://www.uniprot.org/uniprot/{protein_id}.fasta")
        response.raise_for_status()
        return response.text
    except Exception as e:
        st.error(f"UniProt fetch failed: {e}")
        return None

@st.cache_data(ttl=3600)
def fetch_kegg(pathway_id):
    try:
        response = requests.get(f"https://rest.kegg.jp/get/{pathway_id}")
        response.raise_for_status()
        return response.text
    except Exception as e:
        st.error(f"KEGG fetch failed: {e}")
        return None

@st.cache_data(ttl=3600)
def get_snp_info(rsid):
    url = f"https://rest.ensembl.org/variation/human/{rsid}?content-type=application/json"
    response = requests.get(url)
    if response.status_code == 200:
        return response.json()
    else:
        return None

@st.cache_data(ttl=3600)
def get_diseases(rsid):
    url = f"https://www.disgenet.org/api/gda/getDiseases/{rsid}"
    response = requests.get(url)
    if response.status_code == 200:
        diseases = response.json()
        return [disease['diseaseName'] for disease in diseases]
    else:
        return []

@st.cache_data(ttl=3600)
def get_snps_for_gene(gene_name):
    url = f"https://rest.ensembl.org/overlap/region/human/{gene_name}?feature=variation;content-type=application/json"
    response = requests.get(url)
    if response.status_code == 200:
        return response.json()
    else:
        return []

# Main title
st.markdown("<h1 class='title'>🧬SynoBio Toolkit</h1>", unsafe_allow_html=True)

# Horizontal tabs
tabs = st.tabs(["Home", "UniProt Fetch", "BLAST Search", "KEGG Pathway", "PubChem Explorer", "SNP Analysis", "About"])

# Home Tab
with tabs[0]:
    st.markdown("""
    ## Welcome to the SynoBio Toolkit!

    This web server has been developed by a dedicated bioinformatics student to facilitate biological data analysis and research. It integrates multiple tools to help researchers, students, and educators perform various bioinformatics tasks efficiently.

    The server is hosted on a reliable platform to ensure smooth operation and quick access. Designed with user-friendliness in mind, it provides an intuitive interface suitable for users with varying levels of expertise. The platform is regularly updated to include new features and tools that keep up with the latest developments in bioinformatics.

    Feel free to explore the tools available and reach out with feedback or questions. Happy exploring and analyzing!
    """)
    
    st.markdown("""
    <div class="info-box">
        <h3>Why use SynoBio Toolkit?</h3>
        <p>This platform provides researchers, students, and bioinformatics enthusiasts with easy access to powerful biological 
        data analysis tools without requiring programming expertise. All tools are accessible through a simple web interface.</p>
    </div>
    """, unsafe_allow_html=True)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("""
        <div class="info-box">
            <h3>🔬 Bioinformatics Tools</h3>
            <ul>
                <li><strong>BLAST Search</strong>: Compare your sequence against NCBI databases</li>
                <li><strong>UniProt Fetch</strong>: Retrieve protein sequences and annotations</li>
                <li><strong>KEGG Pathway</strong>: Access metabolic pathway information</li>
                <li><strong>PubChem Explorer</strong>: Search chemical compound data</li>
            </ul>
        </div>
        """, unsafe_allow_html=True)
    
    with col2:
        st.markdown("""
        <div class="info-box">
            <h3>🧬 Genetic Analysis Tools</h3>
            <ul>
                <li><strong>SNP Information</strong>: Get detailed data about genetic variants</li>
                <li><strong>Disease Associations</strong>: Explore variant-disease relationships</li>
                <li><strong>Gene Variants</strong>: Find SNPs associated with specific genes</li>
            </ul>
        </div>
        """, unsafe_allow_html=True)
    
    st.markdown("""
    <div class="info-box">
        <h3>🚀 Getting Started</h3>
        <ol>
            <li>Select a tool from the tabs above</li>
            <li>Follow the instructions for each tool</li>
            <li>View and interpret your results</li>
            <li>Download your results for further analysis</li>
        </ol>
    </div>
    """, unsafe_allow_html=True)

    st.markdown("""
    <div class="info-box">
        <h4 style="margin-top:20px;">Key Features:</h4>
        <ul>
            <li><strong>User-friendly interface</strong> requiring no installation</li>
            <li><strong>Direct access</strong> to multiple biological databases</li>
        <li><strong>Quick results</strong> with minimal configuration</li>
        <li><strong>Downloadable results</strong> for further analysis</li>
        <li><strong>Responsive design</strong> works on all devices</li>
    </ul>

    <h4 style="margin-top:20px;">How to Use:</h4>
    <ol>
        <li>Select the tool you need from the navigation tabs</li>
        <li>Enter your query parameters</li>
        <li>Submit your request</li>
        <li>View and interpret the results</li>
        <li>Download results if needed</li>
        Note: All tools are powered by public APIs and databases, providing reliable biological data. 
        Some tools may take several minutes to complete depending on server load and query complexity.
    </ol>
    """, unsafe_allow_html=True)

# BLAST Tab
with tabs[2]:
    st.markdown("<h2 class='tool-header'>🔍 NCBI BLAST Search</h2>", unsafe_allow_html=True)
    
    with st.expander("ℹ️ About this tool"):
        st.markdown("""
        **BLAST (Basic Local Alignment Search Tool)** finds regions of similarity between biological sequences. 
        - Compares your sequence against NCBI's databases
        - Identifies similar sequences from various organisms
        - Helps in gene discovery, protein function prediction, and evolutionary studies
        
        **Output**: XML format results showing sequence matches with alignment details.
        """)
    
    uploaded_file = st.file_uploader("Upload FASTA file", type=["fasta", "fa"])
    blast_program = st.selectbox("Select BLAST program", 
                               ["blastn (nucleotide-nucleotide)", 
                                "blastp (protein-protein)", 
                                "blastx (translated nucleotide-protein)"])
    
    if uploaded_file:
        try:
            record = SeqIO.read(uploaded_file, "fasta")
            st.success(f"Successfully loaded sequence: {record.id}")
            
            if st.button("Run BLAST Search"):
                with st.spinner("Running BLAST (this may take several minutes)..."):
                    program = blast_program.split()[0]
                    blast_results = run_blast(record.seq, program)
                    
                    if blast_results:
                        st.success("BLAST completed successfully!")
                        st.download_button(
                            label="Download BLAST Results (XML)",
                            data=blast_results,
                            file_name="blast_results.xml",
                            help="Download the full BLAST results in XML format"
                        )
        except Exception as e:
            st.error(f"Error reading file: {e}")

# UniProt Tab
with tabs[1]:
    st.markdown("<h2 class='tool-header'>🛡️ UniProt Data Fetcher</h2>", unsafe_allow_html=True)
    
    with st.expander("ℹ️ About this tool"):
        st.markdown("""
        **UniProt** is a comprehensive resource for protein sequence and annotation data.
        - Retrieve protein sequences in FASTA format
        - Access curated protein information
        - Find functional annotations and cross-references
        
        **Output**: FASTA format sequence with header information.
        """)
    
    protein_id = st.text_input("Enter UniProt ID (e.g., P12345):", 
                              placeholder="P12345 or A0A024RBG1",
                              help="Standard UniProt accession numbers only")
    
    if protein_id:
        with st.spinner("Fetching data from UniProt..."):
            fasta_data = fetch_uniprot(protein_id)
            
            if fasta_data:
                st.success("Data retrieved successfully!")
                st.text_area("FASTA Sequence:", fasta_data, height=200)
                
                st.download_button(
                    label="Download FASTA",
                    data=fasta_data,
                    file_name=f"{protein_id}.fasta"
                )
            else:
                st.warning("No data found for this ID. Please check the UniProt ID.")

# KEGG Tab
with tabs[3]:
    st.markdown("<h2 class='tool-header'>🔄 KEGG Pathway Analyzer</h2>", unsafe_allow_html=True)
    
    with st.expander("ℹ️ About this tool"):
        st.markdown("""
        **KEGG** is a database resource for understanding biological pathways.
        - Access metabolic pathways and molecular interactions
        - View pathway maps and associated genes
        - Explore disease pathways and drug targets
        
        **Output**: Textual pathway information and visual pathway map.
        """)
    
    pathway_id = st.text_input(
        "Enter KEGG Pathway ID (e.g., hsa00020):", 
        placeholder="hsa00020 or map00010",
        help="Format: [organism prefix][pathway number]"
    )
    
    # Add a note instructing the user
    st.write("📝 Please ensure the pathway ID is valid and known in KEGG database.")
    
    if pathway_id:
        with st.spinner("Fetching pathway data..."):
            pathway_data = fetch_kegg(pathway_id)
            
            if pathway_data:
                st.success("Pathway data retrieved!")
                st.text_area("Pathway Information:", pathway_data, height=300)
                
                st.image(
                    f"https://www.kegg.jp/kegg/pathway/{pathway_id}.png", 
                    caption=f"Pathway: {pathway_id}",
                    use_container_width=True
                )
            else:
                st.warning("Pathway not found. Check the ID format.")

# PubChem Tab
with tabs[4]:
    st.markdown("<h2 class='tool-header'>⚗️ PubChem Compound Explorer</h2>", unsafe_allow_html=True)

    with st.expander("ℹ️ About this tool"):
        st.markdown("""
        **PubChem** is a chemistry database with information on small molecules.
        - Search chemical compounds by name or SMILES format.
        - View molecular structures and properties.
        - Access chemical identifiers and formulas.
        
        **Output**: Compound information including structure image.
        """)

    # Search by compound name
    compound_name = st.text_input("Enter Compound Name (e.g., Aspirin):", 
                                 placeholder="Aspirin or Caffeine",
                                 help="Common or scientific names of chemical compounds")

    # Search by SMILES string
    smiles_input = st.text_input("Enter SMILES Format:", 
                                placeholder="e.g., CC(=O)OC1=CC=CC=C1C(=O)O", 
                                help="SMILES format for chemical compounds")

    # If compound name is entered
    if compound_name:
        with st.spinner("Searching PubChem..."):
            try:
                compounds = pcp.get_compounds(compound_name, "name")
                
                if compounds:
                    compound = compounds[0]
                    st.success(f"Found: {compound.iupac_name or compound_name}")
                    
                    col1, col2 = st.columns([1, 2])
                    
                    with col1:
                        st.markdown(f"""
                        **CID**: {compound.cid}  
                        **Molecular Formula**: {compound.molecular_formula}  
                        **Molecular Weight**: {compound.molecular_weight} g/mol  
                        **Synonyms**: {', '.join(compound.synonyms[:3])}...
                        """)
                    
                    with col2:
                        st.image(
                            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{compound.cid}/PNG",
                            caption=f"Structure of {compound_name}",
                            width=300
                        )
                else:
                    st.warning("Compound not found. Try alternative names.")
            except Exception as e:
                st.error(f"Search error: {e}")

    # If SMILES string is entered
    if smiles_input:
        with st.spinner("Searching PubChem..."):
            try:
                compounds = pcp.get_compounds(smiles_input, "smiles")
                
                if compounds:
                    compound = compounds[0]
                    st.success(f"Found: {compound.iupac_name or smiles_input}")
                    
                    col1, col2 = st.columns([1, 2])
                    
                    with col1:
                        st.markdown(f"""
                        **CID**: {compound.cid}  
                        **Molecular Formula**: {compound.molecular_formula}  
                        **Molecular Weight**: {compound.molecular_weight} g/mol  
                        **Synonyms**: {', '.join(compound.synonyms[:3])}...
                        """)
                    
                    with col2:
                        st.image(
                            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{compound.cid}/PNG",
                            caption=f"Structure of SMILES: {smiles_input}",
                            width=300
                        )
                else:
                    st.warning("SMILES string not found. Try a valid SMILES.")
            except Exception as e:
                st.error(f"Search error: {e}")


# SNP Analysis Tab
with tabs[5]:
    st.markdown("<h2 class='tool-header'>🧬 SNP Analysis Tools</h2>", unsafe_allow_html=True)
    
    analysis_type = st.radio("Select analysis type:", 
                            ["SNP Information"])
    
    if analysis_type == "SNP Information":
        with st.expander("ℹ️ About this tool"):
            st.markdown("""
            **SNP Information** provides details about genetic variants.
            - Get genomic location and allele information
            - View clinical significance
            - Access evidence sources and synonyms
            
            **Output**: Comprehensive SNP data from Ensembl.
            """)
        
        rsid = st.text_input("Enter SNP ID (e.g., rs429358):", 
                            placeholder="rs429358 or rs7412",
                            help="Standard dbSNP reference SNP IDs (rs numbers)")
        
        if rsid:
            with st.spinner("Fetching SNP data..."):
                snp_data = get_snp_info(rsid)
                
                if snp_data:
                    st.success(f"Data found for {rsid}")
                    
                    subtabs = st.tabs(["Overview", "Clinical Data", "Evidence"])
                    
                    with subtabs[0]:
                        st.markdown(f"""
                        **Name**: {snp_data['name']}  
                        **Location**: {snp_data['mappings'][0]['location']}  
                        **Alleles**: {snp_data['mappings'][0]['allele_string']}  
                        **Ancestral Allele**: {snp_data['mappings'][0]['ancestral_allele']}  
                        **Consequence**: {snp_data['most_severe_consequence']}
                        """)
                    
                    with subtabs[1]:
                        if snp_data['clinical_significance']:
                            st.write("**Clinical Significance:**")
                            for sig in snp_data['clinical_significance']:
                                st.write(f"- {sig}")
                        else:
                            st.write("No clinical significance data available")
                    
                    with subtabs[2]:
                        if snp_data['evidence']:
                            st.write("**Evidence Sources:**")
                            for evidence in snp_data['evidence']:
                                st.write(f"- {evidence}")
                        else:
                            st.write("No evidence sources recorded")
                else:
                    st.warning("SNP not found. Check the ID format.")

# About Tab
with tabs[6]:
    st.markdown("<h2 class='tool-header'>📝 About SynoBio Toolkit</h2>", unsafe_allow_html=True)
    
    st.markdown("""
    ## Overview
    The **SynoBio Toolkit** integrates multiple bioinformatics and genetic analysis tools into a single, user-friendly web interface.
    """)
    
    # Creator Section
    st.markdown("""
    <div class="Creator-section">
        <h3>👩‍🔬 About the Creator</h3>
        <div style="display: flex; align-items: center; gap: 20px;">
            <img src="https://media.licdn.com/dms/image/v2/D4D03AQHGhmKbndo8qg/profile-displayphoto-shrink_400_400/B4DZXhg9kFG4Ag-/0/1743245270040?e=1752105600&v=beta&t=ga1EvqAQtuzFChC4hS21-YA5ZIDgnicXh2P9shgVdqo" 
                 width="120" style="border-radius: 50%; border: 3px solid #6A5ACD;">
            <div>
                <h4>Arati Joshi</h4>
                <p>M.Sc. Bioinformatics Student at DES Pune University</p>
                <p>This web server was developed as a mini project for academic purposes, demonstrating the integration of 
                various bioinformatics tools into a single platform. The Creator is passionate about bioinformatics 
                and computational biology, with interests in genomics and data analysis.</p>
                <p>Arati is dedicated to advancing her skills in biological data interpretation and aims to contribute to research that bridges biology and computer science. She actively participates in bioinformatics communities and workshops to stay updated with the latest trends, aspiring to develop innovative solutions that facilitate scientific discovery.</p>
            </div>
        </div>
    </div>
    """, unsafe_allow_html=True)
    
    # Mentorship Section
    st.markdown("""
    <div class="mentorship-section">
        <h3>👨‍🏫 Mentorship</h3>
        <div style="display: flex; align-items: center; gap: 20px;">
            <img src="https://media.licdn.com/dms/image/v2/D5603AQF9gsU7YBjWVg/profile-displayphoto-shrink_400_400/B56ZZI.WrdH0Ag-/0/1744981029051?e=1752105600&v=beta&t=F4QBDSEgjUvnBS00xPkKqPTLI0jQaMpYefaOzARY1Yg" 
                width="120" style="border-radius: 50%; border: 3px solid #6A5ACD;">
            <div>
                <h4>Dr. Kushagra Kashyap</h4>
                <p>Assistant Professor (Bioinformatics), Department of Life Sciences, School of Science and Mathematics, DES Pune University</p>
                <p>This project was developed under the guidance of Dr. Kashyap, who provided valuable insights and mentorship 
                throughout the development process. His expertise in bioinformatics and computational biology was instrumental 
                in shaping this project.</p>
                <a href="https://www.linkedin.com/in/dr-kushagra-kashyap-b230a3bb" target="_blank">🔗 Connect on LinkedIn</a>
            </div>
        </div>
    </div>
    """, unsafe_allow_html=True)

    
    # Web Server Info Section
    st.markdown("""
    # About the Web Server

    **The SynoBio Toolkit** is designed to provide researchers and students with easy access to essential bioinformatics tools without requiring programming knowledge. The server integrates multiple public APIs and databases to deliver comprehensive biological data analysis capabilities.
    """)

    with st.expander("🔧 Tools Included"):
        st.markdown("""
        - **UniProt Fetch**: Protein sequence retrieval
        - **BLAST Search**: Sequence similarity searching
        - **KEGG Pathway**: Metabolic pathway analysis
        - **PubChem Explorer**: Chemical compound information
        - **SNP Analysis**: Genetic variant exploration
        """)

    with st.expander("📚 Data Sources"):
        st.markdown("""
        - UniProt (Universal Protein Resource)
        - NCBI BLAST (National Center for Biotechnology Information)
        - KEGG (Kyoto Encyclopedia of Genes and Genomes)
        - PubChem (NIH Chemical Database)
        - Ensembl (Genome Database)
        - DisGeNET (Disease-Gene-Variant Network)
        """)

    st.markdown("""
    ### Technical Details
    - Python 3
    - Streamlit (Web Framework)
    - Biopython (Bioinformatics Library)
    - PubChemPy (Chemical Data Access)
    - Requests (API Communication)

    ---

    ### Feedback & Contact

    **Connect with me**  
    Email: [joshiarati810@gmail.com](mailto:joshiarati810@gmail.com)  
    LinkedIn: [www.linkedin.com/in/arati-joshi-a2719819a](https://www.linkedin.com/in/arati-joshi-a2719819a)  
    GitHub: [github.com/AratiJ17](https://github.com/AratiJ17)

    <br>
    Thank you for using the SynoBio Toolkit!

    <div style="margin-top:30px; font-size: 0.9em; color: gray;" class="footer">
        2023 SynoBio Toolkit | DES Pune University | Academic Project
    </div>
    """, unsafe_allow_html=True)
