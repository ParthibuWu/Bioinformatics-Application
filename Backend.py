from fastapi import FastAPI, File, UploadFile, HTTPException
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import io
import gzip
import bz2
from typing import List, Dict, Set, Optional

app = FastAPI(title="Bioinformatics Processing API")


def calculate_gc_content(seq, mode="raw") -> float:
    """Calculate GC content"""
    seq = seq.upper()
    gc = 0
    canonical = 0

    for base in seq:
        if base in ("G", "C"):
            gc += 1
            canonical += 1
        elif base in ("A", "T"):
            canonical += 1

    if mode == "raw":
        denominator = len(seq)
    elif mode == "canonical":
        denominator = canonical
    else:
        raise ValueError("Mode must be either 'raw' or 'canonical'")

    if denominator == 0:
        return 0.0

    return round((gc / denominator) * 100, 2)


def detect_compression(filename: str) -> str:
    """Detect compression type from filename"""
    filename_lower = filename.lower()

    if filename_lower.endswith(('.gz', '.gzip')):
        return 'gzip'
    elif filename_lower.endswith(('.bz2', '.bzip2')):
        return 'bzip2'
    else:
        return 'none'


def detect_file_format(filename: str) -> str:
    """Detect biological file format from extension"""
    name = filename.lower()

    # Remove compression extensions
    for ext in ['.gz', '.gzip', '.bz2', '.bzip2']:
        if name.endswith(ext):
            name = name[:-len(ext)]

    # Check biological format
    if name.endswith(('.fa', '.fasta', '.fna')):
        return 'fasta'
    elif name.endswith(('.fq', '.fastq')):
        return 'fastq'
    elif name.endswith(('.gb', '.gbk', '.genbank')):
        return 'genbank'
    elif name.endswith('.embl'):
        return 'embl'
    else:
        raise ValueError(f"Unsupported file format: {filename}")


def get_text_handle(content: bytes, compression: str):
    """
    Convert bytes to text-mode file handle
    Returns a file-like object that SeqIO can read
    """
    if compression == 'gzip':
        binary_io = io.BytesIO(content)
        return gzip.open(binary_io, 'rt')  # 'rt' = read text mode

    elif compression == 'bzip2':
        binary_io = io.BytesIO(content)
        return bz2.open(binary_io, 'rt')

    else:
        # Uncompressed - decode to string
        text_content = content.decode('utf-8')
        return io.StringIO(text_content)


def process_fasta_stream(file_handle, file_format: str = "fasta") -> List[Dict]:
    """Process FASTA/GenBank files"""
    records_data = []

    for record in SeqIO.parse(file_handle, file_format):
        record_info = {
            "ID": record.id,
            "Description": record.description,
            "Sequence": str(record.seq),
            "Length": len(record),
            "GC_content": calculate_gc_content(record.seq),
            "Last_base": str(record.seq[-1]) if len(record.seq) > 0 else "",
            "First_base": str(record.seq[0]) if len(record.seq) > 0 else ""
        }
        records_data.append(record_info)

    return records_data


def process_fastq_stream(file_handle, wanted_ids: Optional[Set[str]] = None) -> List[Dict]:
    """Process FASTQ files with optional filtering"""
    if wanted_ids is None:
        wanted_ids = set()

    records_data = []

    # Read content for FastqGeneralIterator
    content = file_handle.read()
    if isinstance(content, bytes):
        content = content.decode('utf-8')

    for title, seq, qual in FastqGeneralIterator(io.StringIO(content)):
        record_id = title.split(None, 1)[0]

        if wanted_ids and record_id not in wanted_ids:
            continue

        record_info = {
            "ID": record_id,
            "Title": title,
            "Sequence": seq,
            "Quality": qual,
            "Length": len(seq),
            "GC_content": calculate_gc_content(seq),
            "Avg_quality": sum(ord(q) - 33 for q in qual) / len(qual) if qual else 0
        }
        records_data.append(record_info)

    return records_data


@app.post("/sequences/process/")
async def process_sequences_universal(file: UploadFile = File(...)):
    """
    Universal endpoint - handles FASTA, FASTQ, GenBank
    Supports compressed (.gz, .bz2) files
    """
    try:
        # Detect compression and format
        compression = detect_compression(file.filename)
        file_format = detect_file_format(file.filename)

        # CRITICAL FIX: Read file content first (await the async operation)
        content = await file.read()

        # Then get text handle from bytes (synchronous operation)
        handle = get_text_handle(content, compression)

        # Process based on format
        if file_format == 'fastq':
            sequences = process_fastq_stream(handle)
        else:
            sequences = process_fasta_stream(handle, file_format=file_format)

        return {
            "filename": file.filename,
            "format": file_format,
            "compression": compression,
            "total_sequences": len(sequences),
            "total_bases": sum(s["Length"] for s in sequences),
            "sequences": sequences
        }

    except ValueError as e:
        raise HTTPException(400, str(e))
    except Exception as e:
        raise HTTPException(500, f"Error processing file: {str(e)}")


@app.post("/sequences/stats/")
async def get_sequence_stats(file: UploadFile = File(...)):
    """Get statistics only - more memory efficient"""
    try:
        compression = detect_compression(file.filename)
        file_format = detect_file_format(file.filename)

        # CRITICAL FIX: await file.read() first
        content = await file.read()
        handle = get_text_handle(content, compression)

        # Calculate stats
        total_length = 0
        total_gc = 0
        count = 0

        for record in SeqIO.parse(handle, file_format):
            total_length += len(record)
            gc_content = calculate_gc_content(record.seq)
            total_gc += gc_content
            count += 1

        return {
            "filename": file.filename,
            "format": file_format,
            "compression": compression,
            "total_sequences": count,
            "total_bases": total_length,
            "average_length": round(total_length / count, 2) if count > 0 else 0,
            "average_gc_content": round(total_gc / count, 2) if count > 0 else 0
        }

    except ValueError as e:
        raise HTTPException(400, str(e))
    except Exception as e:
        raise HTTPException(500, f"Error processing file: {str(e)}")


@app.post("/fastq/filter/")
async def filter_fastq(
        file: UploadFile = File(...),
        filter_ids: Optional[str] = None
):
    """
    FASTQ-specific endpoint with ID filtering
    Supports compressed files
    """
    try:
        compression = detect_compression(file.filename)
        file_format = detect_file_format(file.filename)

        if file_format != 'fastq':
            raise HTTPException(400, "This endpoint only accepts FASTQ files")

        # Parse filter IDs
        wanted_ids = None
        if filter_ids:
            wanted_ids = set(id.strip() for id in filter_ids.split(","))

        # CRITICAL FIX: await file.read() first
        content = await file.read()
        handle = get_text_handle(content, compression)

        sequences = process_fastq_stream(handle, wanted_ids=wanted_ids)

        return {
            "filename": file.filename,
            "compression": compression,
            "total_sequences": len(sequences),
            "filtered": bool(wanted_ids),
            "filter_count": len(wanted_ids) if wanted_ids else 0,
            "sequences": sequences
        }

    except ValueError as e:
        raise HTTPException(400, str(e))
    except Exception as e:
        raise HTTPException(500, f"Error processing file: {str(e)}")


@app.get("/")
async def root():
    return {
        "message": "Bioinformatics Processing API",
        "supported_formats": ["FASTA (.fa, .fasta)", "FASTQ (.fq, .fastq)", "GenBank (.gb, .gbk)", "EMBL (.embl)"],
        "supported_compression": ["gzip (.gz)", "bzip2 (.bz2)", "uncompressed"],
        "endpoints": {
            "POST /sequences/process/": "Process any sequence file",
            "POST /sequences/stats/": "Get statistics only (memory efficient)",
            "POST /fastq/filter/": "Filter FASTQ by sequence IDs"
        }
    }