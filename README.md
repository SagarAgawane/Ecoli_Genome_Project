# 🧬 E. coli B REL606 Genome Assembly & Annotation  
### **Genome Assembly, Polishing, and Annotation using Illumina Data**
[![GitHub](https://badgen.net/badge/icon/github?icon=github&label)](https://github.com/SagarAgawane/Ecoli_Genome_Project)

## 📌 **Project Overview**  
This project focuses on the **assembly, polishing, and annotation** of the *Escherichia coli* B strain **REL606** genome using **Illumina sequencing reads from NCBI SRA (SRR2584863)**.  
**Objective:** Generate a high-quality, annotated genome using bioinformatics tools.  

---

## 📂 **Project Structure**

- **annotation_output_v4/** → *Final genome annotation (Prokka output)*
- **ecoli_assembly/** → *Genome assembly (SPAdes output)*
- **ecoli_pilon/** → *Polished genome (Pilon output)*
- **corrected_ecoli.fasta** → *Final genome sequence*
- **corrected_ecoli.sqn** → *Final NCBI submission file (if needed)*
- **corrected_ecoli.val** → *Validation report*
- **template.sbt** → *Submission metadata (NCBI-style)*
- **biosample_result.xml** → *Metadata from NCBI BioSample*
- **README.md** → *Project documentation*

---

## 📌 Workflow  
### **1️⃣ Data Download & Quality Control**  
- Raw sequencing data was downloaded from **NCBI SRA** (SRR2584863).  
- Quality control and **adapter trimming** were performed using **FastQC** and **Trimmomatic**.
### **2️⃣ Genome Assembly & Polishing**  
- The genome was assembled using **SPAdes**.  
- Sorting and indexing of BAM file was performed using **SAMtools**, and the assembly was polished using **Pilon**,which internally aligned the reads.

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

⚠ **Important Note:**  
This project was conducted **for educational purposes only** and **has not been submitted to NCBI GenBank** because it is based on publicly available sequencing data from **NCBI SRA (SRR2584863)**.  

The `.sqn` file was generated using `tbl2asn` as a learning exercise, but **this genome is not officially deposited in any public database**.  

### ❌ **Why Wasn't This Uploaded to NCBI?**
- This project **uses sequencing data from NCBI SRA**, which is already public.  
- Submitting this assembly would be **redundant and against NCBI policies**.  
- The genome annotation was performed **for learning purposes only**, not for official publication.  

If you wish to submit a genome to NCBI, ensure that:  
- You **own the sequencing data** or have explicit permission to use it.  
- The metadata and annotations **comply with NCBI's submission guidelines**.  
- The genome assembly **meets NCBI quality standards**.  

For official genome submissions, refer to the [NCBI Submission Guide](https://www.ncbi.nlm.nih.gov/genbank/submit/).  

---

## 🔍 Visualization  
- The annotated genome can be viewed using:  
  - **Artemis** (for `.gbk` files)  
  - **IGV** (for `.gff` + `.fasta`)  
  - **NCBI Genome Workbench** (for `.gbk` online)  

```bash
# Example: Open annotation in Artemis
artemis annotation_output_v4/ecoli_v4.gbk
```
---

## 👨‍💻 Author & Contact  
**Name:** Sagar Agavane  
**GitHub:** [SagarAgawane](https://github.com/SagarAgawane)  
**Email:** sagar0203agavane@gmail.com  

