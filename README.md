# 🧬 E. coli B REL606 Genome Assembly & Annotation  
### **Genome Assembly, Polishing, and Annotation using Illumina Data**
[![GitHub](https://badgen.net/badge/icon/github?icon=github&label)](https://github.com/SagarAgawane/Ecoli_Genome_Project)

## 📌 **Project Overview**  
This project focuses on the **assembly, polishing, and annotation** of the *Escherichia coli* B strain **REL606** genome using **Illumina sequencing reads from NCBI SRA (SRR2584863)**.  
**Objective:** Generate a high-quality, annotated genome using bioinformatics tools.  

---

## 📂 **Project Structure**
Ecoli_Genome_Project/ │── annotation_output_v4/ # Final genome annotation (Prokka output) │── ecoli_assembly/ # Genome assembly (SPAdes output) │── ecoli_pilon/ # Polished genome (Pilon output) │── corrected_ecoli.fasta # Final genome sequence │── corrected_ecoli.sqn # Final NCBI submission file │── corrected_ecoli.val # Validation report │── template.sbt # Submission metadata (NCBI) │── biosample_result.xml # Metadata from NCBI BioSample │── README.md # Project documentation

---

## 📌 Workflow  
### **1️⃣ Data Download & Quality Control**  
- Raw sequencing data was downloaded from **NCBI SRA** (SRR2584863).  
- Quality control and **adapter trimming** were performed using **FastQC** and **Trimmomatic**.

### **2️⃣ Genome Assembly & Polishing**  
- The genome was assembled using **SPAdes**.  
- Read alignment was performed using **BWA**, and the assembly was polished using **Pilon**.

### **3️⃣ Genome Annotation**  
- The polished genome was annotated using **Prokka**.  
- Functional features (CDS, rRNAs, tRNAs) were identified and verified.

### **4️⃣ NCBI Submission Preparation**  
- The final genome sequence was **validated using `tbl2asn`**.  
- The **`.sqn` file** was generated for submission to **NCBI GenBank**.

---

## 📜 Tools & Commands Used  
| **Tool**  | **Purpose** | **Command Example** |
|-----------|------------|--------------------|
| `FastQC` | Quality Check | `fastqc SRR2584863_1.fastq.gz` |
| `Trimmomatic` | Read Trimming | `trimmomatic PE -threads 4 ...` |
| `SPAdes` | Genome Assembly | `spades.py --careful -1 trimmed_1P.fastq.gz -2 trimmed_2P.fastq.gz -o ecoli_assembly` |
| `BWA` | Read Mapping | `bwa mem ecoli_assembly/scaffolds.fasta trimmed_1P.fastq.gz trimmed_2P.fastq.gz > aligned_reads.sam` |
| `SAMtools` | BAM Processing | `samtools view -bS aligned_reads.sam | samtools sort -o aligned_reads.sorted.bam` |
| `Pilon` | Genome Polishing | `pilon --genome ecoli_assembly/scaffolds.fasta --bam aligned_reads.sorted.bam --outdir ecoli_pilon --fix all --threads 4` |
| `Prokka` | Genome Annotation | `prokka --outdir annotation_output_v4 --prefix ecoli_v4 corrected_ecoli.fasta` |
| `tbl2asn` | NCBI Submission Validation | `tbl2asn -i corrected_ecoli.fasta -s template.sbt -M n -J` |

---

## 🗃️ Final Genome Statistics  
| Feature | Count |
|---------|-------|
| **Contigs** | 107 |
| **Genome Size** | 4,563,017 bp |
| **GC Content** | 50.79% |
| **CDS (Protein-Coding Genes)** | 4,261 |
| **tRNAs** | 72 |
| **rRNAs** | 5 |
| **tmRNA** | 1 |

---

## 📤 Submission to NCBI  
The final genome annotation was validated and is **ready for submission** to NCBI GenBank.  
✅ **Submission File:** [`corrected_ecoli.sqn`](corrected_ecoli.sqn)  
✅ **Validation Report:** [`corrected_ecoli.val`](corrected_ecoli.val)  

Submission can be done via the **[NCBI Submission Portal](https://submit.ncbi.nlm.nih.gov/genbank/)**.

---

## 🔍 Visualization  
- The annotated genome can be viewed using:  
  - **Artemis** (for `.gbk` files)  
  - **IGV** (for `.gff` + `.fasta`)  
  - **NCBI Genome Workbench** (for `.gbk` online)  

```bash
# Example: Open annotation in Artemis
artemis annotation_output_v4/ecoli_v4.gbk

👨‍💻 Author & Contact
Name: Sagar Agawane
GitHub: Sagar Agavane
Email: sagar0203agavane@gmail.com
