import streamlit as st
import requests
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from io import StringIO
import json

# Configure page
st.set_page_config(
    page_title="Bioinformatics Processing Tool",
    page_icon="🧬",
    layout="wide"
)

# API base URL - change if your API runs on different port
API_BASE_URL = "http://localhost:8000"

# Custom CSS for better styling
st.markdown("""
    <style>
    .stAlert {
        padding: 1rem;
        margin: 1rem 0;
    }
    .metric-card {
        background-color: #f0f2f6;
        padding: 1rem;
        border-radius: 0.5rem;
        margin: 0.5rem 0;
    }
    </style>
""", unsafe_allow_html=True)


def check_api_health():
    """Check if API is running"""
    try:
        response = requests.get(f"{API_BASE_URL}/")
        return response.status_code == 200
    except:
        return False


def process_sequences(file, endpoint="/sequences/process/"):
    """Send file to API and get results"""
    try:
        files = {"file": (file.name, file.getvalue(), file.type)}
        response = requests.post(f"{API_BASE_URL}{endpoint}", files=files)

        if response.status_code == 200:
            return response.json(), None
        else:
            return None, f"Error {response.status_code}: {response.text}"
    except Exception as e:
        return None, f"Connection error: {str(e)}"


def process_with_filter(file, filter_ids):
    """Process FASTQ with ID filtering"""
    try:
        files = {"file": (file.name, file.getvalue(), file.type)}
        params = {"filter_ids": filter_ids} if filter_ids else {}

        response = requests.post(
            f"{API_BASE_URL}/fastq/filter/",
            files=files,
            params=params
        )

        if response.status_code == 200:
            return response.json(), None
        else:
            return None, f"Error {response.status_code}: {response.text}"
    except Exception as e:
        return None, f"Connection error: {str(e)}"


def get_stats_only(file):
    """Get statistics without full sequence data"""
    try:
        files = {"file": (file.name, file.getvalue(), file.type)}
        response = requests.post(f"{API_BASE_URL}/sequences/stats/", files=files)

        if response.status_code == 200:
            return response.json(), None
        else:
            return None, f"Error {response.status_code}: {response.text}"
    except Exception as e:
        return None, f"Connection error: {str(e)}"


def plot_gc_distribution(sequences):
    """Create GC content distribution plot"""
    gc_values = [seq["GC_content"] for seq in sequences]

    fig = px.histogram(
        x=gc_values,
        nbins=30,
        title="GC Content Distribution",
        labels={"x": "GC Content (%)", "y": "Count"},
        color_discrete_sequence=["#636EFA"]
    )

    fig.update_layout(
        showlegend=False,
        height=400
    )

    return fig


def plot_length_distribution(sequences):
    """Create sequence length distribution plot"""
    lengths = [seq["Length"] for seq in sequences]

    fig = px.histogram(
        x=lengths,
        nbins=30,
        title="Sequence Length Distribution",
        labels={"x": "Length (bp)", "y": "Count"},
        color_discrete_sequence=["#EF553B"]
    )

    fig.update_layout(
        showlegend=False,
        height=400
    )

    return fig


def plot_base_composition(sequences):
    """Plot first and last base composition"""
    first_bases = [seq["First_base"] for seq in sequences if "First_base" in seq]
    last_bases = [seq["Last_base"] for seq in sequences if "Last_base" in seq]

    # Count bases
    from collections import Counter
    first_counts = Counter(first_bases)
    last_counts = Counter(last_bases)

    fig = go.Figure(data=[
        go.Bar(name='First Base', x=list(first_counts.keys()), y=list(first_counts.values())),
        go.Bar(name='Last Base', x=list(last_counts.keys()), y=list(last_counts.values()))
    ])

    fig.update_layout(
        title="First and Last Base Composition",
        xaxis_title="Base",
        yaxis_title="Count",
        barmode='group',
        height=400
    )

    return fig


def main():
    # Header
    st.title("🧬 Bioinformatics Sequence Processing Tool")
    st.markdown("---")

    # Check API health
    if not check_api_health():
        st.error("⚠️ Cannot connect to API. Make sure the FastAPI server is running on http://localhost:8000")
        st.info("Run: `uvicorn Ship:app --reload` or `python -m uvicorn Ship:app --reload`")
        return

    st.success("✅ Connected to API")

    # Sidebar for configuration
    with st.sidebar:
        st.header("⚙️ Configuration")

        processing_mode = st.radio(
            "Processing Mode",
            ["Full Processing", "Statistics Only", "FASTQ Filtering"]
        )

        st.markdown("---")
        st.subheader("Supported Formats")
        st.markdown("""
        **Sequence Formats:**
        - FASTA (.fa, .fasta, .fna)
        - FASTQ (.fq, .fastq)
        - GenBank (.gb, .gbk)
        - EMBL (.embl)

        **Compression:**
        - gzip (.gz)
        - bzip2 (.bz2)
        - uncompressed
        """)

    # Main content area
    uploaded_file = st.file_uploader(
        "Upload your sequence file",
        type=["fa", "fasta", "fna", "fq", "fastq", "gb", "gbk", "genbank", "embl", "gz", "bz2"],
        help="Supports FASTA, FASTQ, GenBank, EMBL formats (compressed or uncompressed)"
    )

    if uploaded_file is not None:
        # Show file info
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Filename", uploaded_file.name)
        with col2:
            st.metric("File Size", f"{uploaded_file.size / 1024:.2f} KB")
        with col3:
            st.metric("Type", uploaded_file.type or "Unknown")

        st.markdown("---")

        # Processing based on mode
        if processing_mode == "Statistics Only":
            if st.button("📊 Calculate Statistics", type="primary"):
                with st.spinner("Processing..."):
                    result, error = get_stats_only(uploaded_file)

                if error:
                    st.error(error)
                else:
                    st.success("✅ Statistics calculated successfully!")

                    # Display stats
                    col1, col2, col3, col4 = st.columns(4)

                    with col1:
                        st.metric("Total Sequences", result["total_sequences"])
                    with col2:
                        st.metric("Total Bases", f"{result['total_bases']:,}")
                    with col3:
                        st.metric("Avg Length", f"{result['average_length']:.2f} bp")
                    with col4:
                        st.metric("Avg GC%", f"{result['average_gc_content']:.2f}%")

                    # Show file details
                    st.info(f"📁 Format: {result['format'].upper()} | Compression: {result['compression'].upper()}")

        elif processing_mode == "FASTQ Filtering":
            filter_ids = st.text_area(
                "Filter by Sequence IDs (comma-separated)",
                placeholder="seq1, seq2, seq3",
                help="Leave empty to process all sequences"
            )

            if st.button("🔍 Process & Filter", type="primary"):
                with st.spinner("Processing..."):
                    result, error = process_with_filter(uploaded_file, filter_ids)

                if error:
                    st.error(error)
                else:
                    st.success("✅ Processing complete!")

                    # Display summary
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        st.metric("Total Sequences", result["total_sequences"])
                    with col2:
                        st.metric("Filtered", "Yes" if result["filtered"] else "No")
                    with col3:
                        st.metric("Filter Count", result.get("filter_count", 0))

                    # Display sequences
                    if result["sequences"]:
                        st.subheader("📋 Filtered Sequences")

                        # Convert to DataFrame
                        df = pd.DataFrame(result["sequences"])

                        # Show interactive table
                        st.dataframe(
                            df,
                            use_container_width=True,
                            height=400
                        )

                        # Download button
                        csv = df.to_csv(index=False)
                        st.download_button(
                            label="⬇️ Download Results (CSV)",
                            data=csv,
                            file_name="filtered_sequences.csv",
                            mime="text/csv"
                        )

                        # Visualizations
                        st.subheader("📊 Visualizations")

                        col1, col2 = st.columns(2)

                        with col1:
                            st.plotly_chart(plot_gc_distribution(result["sequences"]), use_container_width=True)

                        with col2:
                            st.plotly_chart(plot_length_distribution(result["sequences"]), use_container_width=True)

                        if "Avg_quality" in df.columns:
                            st.subheader("Quality Scores Distribution")
                            fig = px.histogram(
                                df,
                                x="Avg_quality",
                                nbins=30,
                                title="Average Quality Score Distribution"
                            )
                            st.plotly_chart(fig, use_container_width=True)

        else:  # Full Processing
            if st.button("🚀 Process File", type="primary"):
                with st.spinner("Processing sequences..."):
                    result, error = process_sequences(uploaded_file)

                if error:
                    st.error(error)
                else:
                    st.success("✅ Processing complete!")

                    # Display summary metrics
                    col1, col2, col3, col4 = st.columns(4)

                    with col1:
                        st.metric("Format", result["format"].upper())
                    with col2:
                        st.metric("Compression", result["compression"].upper())
                    with col3:
                        st.metric("Sequences", result["total_sequences"])
                    with col4:
                        st.metric("Total Bases", f"{result['total_bases']:,}")

                    # Tabs for different views
                    tab1, tab2, tab3, tab4 = st.tabs(["📋 Data Table", "📊 Visualizations", "🔍 Details", "📥 Export"])

                    with tab1:
                        # Convert to DataFrame
                        df = pd.DataFrame(result["sequences"])

                        # Show summary stats
                        st.subheader("Summary Statistics")
                        col1, col2, col3 = st.columns(3)
                        with col1:
                            st.metric("Mean Length", f"{df['Length'].mean():.2f} bp")
                        with col2:
                            st.metric("Mean GC%", f"{df['GC_content'].mean():.2f}%")
                        with col3:
                            st.metric("Std Dev GC%", f"{df['GC_content'].std():.2f}%")

                        # Interactive table
                        st.dataframe(
                            df,
                            use_container_width=True,
                            height=500
                        )

                    with tab2:
                        st.subheader("Sequence Statistics")

                        col1, col2 = st.columns(2)

                        with col1:
                            st.plotly_chart(plot_gc_distribution(result["sequences"]), use_container_width=True)

                        with col2:
                            st.plotly_chart(plot_length_distribution(result["sequences"]), use_container_width=True)

                        # Base composition if available
                        if "First_base" in result["sequences"][0]:
                            st.plotly_chart(plot_base_composition(result["sequences"]), use_container_width=True)

                    with tab3:
                        st.subheader("Detailed Sequence Information")

                        # Sequence selector
                        sequence_index = st.selectbox(
                            "Select sequence to view",
                            range(len(result["sequences"])),
                            format_func=lambda x: result["sequences"][x]["ID"]
                        )

                        selected_seq = result["sequences"][sequence_index]

                        # Display details
                        col1, col2 = st.columns(2)

                        with col1:
                            st.markdown(f"**ID:** {selected_seq['ID']}")
                            st.markdown(f"**Description:** {selected_seq.get('Description', 'N/A')}")
                            st.markdown(f"**Length:** {selected_seq['Length']} bp")
                            st.markdown(f"**GC Content:** {selected_seq['GC_content']}%")

                            if "First_base" in selected_seq:
                                st.markdown(f"**First Base:** {selected_seq['First_base']}")
                                st.markdown(f"**Last Base:** {selected_seq['Last_base']}")

                        with col2:
                            # Show sequence (truncated if too long)
                            seq = selected_seq["Sequence"]
                            if len(seq) > 500:
                                st.text_area(
                                    "Sequence (first 500 bp)",
                                    seq[:500] + "...",
                                    height=200
                                )
                                st.caption(f"Full sequence length: {len(seq)} bp")
                            else:
                                st.text_area(
                                    "Sequence",
                                    seq,
                                    height=200
                                )

                    with tab4:
                        st.subheader("Export Options")

                        # CSV export
                        csv = df.to_csv(index=False)
                        st.download_button(
                            label="⬇️ Download as CSV",
                            data=csv,
                            file_name=f"{uploaded_file.name}_processed.csv",
                            mime="text/csv"
                        )

                        # JSON export
                        json_str = json.dumps(result, indent=2)
                        st.download_button(
                            label="⬇️ Download as JSON",
                            data=json_str,
                            file_name=f"{uploaded_file.name}_processed.json",
                            mime="application/json"
                        )

                        # Excel export option
                        st.info("💡 Tip: Import the CSV into Excel for further analysis")


if __name__ == "__main__":
    main()