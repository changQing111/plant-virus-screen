# Workflow diagram

The arrows show evidence and data dependencies, not concurrent execution. Read cleaning and transcriptome assembly occur upstream. All candidates are retained, including cellular-best conflicts.

```mermaid
flowchart TD
    A["Assembly FASTA + optional reads"]
    B["Discover & merge candidates<br/>BLASTN · DIAMOND · RdRp HMM"]
    C["Evaluate evidence<br/>Competitive search + optional read mapping"]
    D["Annotate & grade<br/>ICTV VMR + evidence grades A–D"]
    E["Validate & report<br/>Candidate tables · summaries · plots"]
    R["Existing enhanced results"]
    A --> B --> C --> D --> E
    R -.-> D
    classDef input fill:#e8f1fb,stroke:#3974a8,color:#173652;
    classDef process fill:#e8f5ef,stroke:#378268,color:#173e30;
    classDef result fill:#eee9fa,stroke:#8063b2,color:#3e2d60;
    class A,R input;
    class B,C,D process;
    class E result;
```

Detailed thresholds and grade rules are documented in the README. D takes precedence over A/B/C. The evidence categories are not statistically independent. Plot selection never removes candidates from the report. Existing-result mode bypasses discovery and mapping; it uses the evidence already present in the completed source run. HMM profiles are selected by RNA-polymerase keywords and do not comprehensively cover all virus groups.
