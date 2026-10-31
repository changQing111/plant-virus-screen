# Workflow diagram

The arrows show evidence and data dependencies, not concurrent execution. Read cleaning and transcriptome assembly occur upstream. All candidates are retained, including cellular-best conflicts.

```mermaid
flowchart TD
    AS["Assembled transcriptome FASTA<br/>Full assembly; unique contig IDs"]
    RD["Optional paired-end reads"]
    DB["Viral nucleotide and protein databases<br/>Pfam RNA polymerase profiles"]
    IN["1. Validate inputs and resources<br/>Initialize manifests and command logs"]
    NT["2. BLASTN<br/>Length >= 300 nt; E <= 1e-10<br/>Bit score >= 200"]
    AA["3. DIAMOND BLASTX<br/>Very sensitive; length >= 100 aa<br/>E <= 1e-5; bit score >= 80"]
    HM["4. Six-frame translation + HMMER<br/>Peptides >= 80 aa; domain i-E <= 1e-5<br/>Model coverage >= 0.25"]
    UNION["5. Union of nt / protein / RdRP evidence<br/>Retain candidates; audit every contig"]
    CP["Competitive protein database<br/>Exact-ID taxonomy mapping"]
    COMP["6. Competitive DIAMOND BLASTX<br/>Classify highest-bit-score hit<br/>viral / cellular / unclassified"]
    MAPQ{"7. Paired reads supplied<br/>and candidates present?"}
    MAP["BWA-MEM against full assembly<br/>samtools BAM checks and coverage<br/>All alignments and MAPQ >= 20"]
    SKIP["Mapping skipped<br/>Missing coverage is not negative evidence"]
    REF["Reference FASTAs + ICTV VMR workbook"]
    ANN["8. Reference and ICTV VMR annotation<br/>Reference-related names, genome and host fields"]
    OLD["Completed enhanced results<br/>--existing-results<br/>Validated source files; optional coverage"]
    DQ{"9. Competitive best hit cellular?"}
    D["D: Conflicting"]
    TWO{"At least two evidence categories?"}
    C["C: Provisional"]
    AQ{"Viral competitive best hit AND<br/>MAPQ >= 20 breadth >= 0.80 AND<br/>mean depth >= 5?"}
    A["A: Supported"]
    B["B: Supported; incomplete validation"]
    SUM["Reference-related group summaries<br/>Keep all candidate IDs, grades and reasons"]
    PLOT["10. Alignment HTML / SVG plots<br/>Default: A/B with applicable alignments<br/>Optional PDF export"]
    VAL["Validate IDs, order, counts and group totals<br/>Record hashes, thresholds and run status"]
    OUT["Candidate table + group summaries<br/>Plots + report + validation.json"]
    NOTE["Heuristic review priorities, not diagnoses<br/>Reference host does not establish plant infection<br/>Groups are not confirmed species counts"]
    AS --> IN
    DB --> IN
    IN --> NT
    IN --> AA
    IN --> HM
    NT --> UNION
    AA --> UNION
    HM --> UNION
    UNION --> COMP
    CP --> COMP
    COMP --> MAPQ
    RD -.-> MAPQ
    MAPQ -->|Yes| MAP
    AS -.-> MAP
    MAPQ -->|No| SKIP
    MAP --> ANN
    SKIP --> ANN
    REF --> ANN
    OLD --> ANN
    ANN --> DQ
    DQ -->|Yes| D
    DQ -->|No| TWO
    TWO -->|No| C
    TWO -->|Yes| AQ
    AQ -->|Yes| A
    AQ -->|No| B
    A --> SUM
    B --> SUM
    C --> SUM
    D --> SUM
    SUM --> PLOT
    PLOT --> VAL
    VAL --> OUT
    OUT -.-> NOTE
    classDef input fill:#e8f1fb,stroke:#3974a8,color:#173652;
    classDef process fill:#e8f5ef,stroke:#378268,color:#173e30;
    classDef decision fill:#fff4d6,stroke:#b58a25,color:#57400d;
    classDef conflict fill:#fbe8e8,stroke:#b85858,color:#622525;
    classDef result fill:#eee9fa,stroke:#8063b2,color:#3e2d60;
    class AS,RD,DB,CP,REF,OLD input;
    class IN,NT,AA,HM,UNION,COMP,MAP,SKIP,ANN,SUM,PLOT,VAL process;
    class MAPQ,DQ,TWO,AQ decision;
    class D,NOTE conflict;
    class A,B,C,OUT result;
```

Default grade thresholds are shown. D takes precedence over A/B/C. The evidence categories are not statistically independent. Plot selection never removes candidates from the report. Existing-result mode bypasses discovery and mapping; it uses the evidence already present in the completed source run. HMM profiles are selected by RNA-polymerase keywords and do not comprehensively cover all virus groups.
