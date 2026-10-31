# plant-virus-screen

以 enhanced 为主干的植物病毒候选发现和证据汇总流程（开发版 0.1.0）。

入口为 screen_virus_res.sh；Python 入口为 unified_pipeline.py。原始流程和数据未修改。

## 运行环境
Python 3（本机以 /data/changq/miniforge3/bin/python 验证）、Biopython。新发现阶段还需 BLASTN、DIAMOND、HMMER，以及需要回贴时的 BWA/samtools。HTML/SVG 绘图使用已有模块；PDF 可选，依赖转换工具，尚未验证。

## 已有结果后处理
```bash
PYTHON_BIN=/data/changq/miniforge3/bin/python bash screen_virus_res.sh \
  --existing-results ../LGTZ_enhanced_res --out new_report \
  --vmr /path/to/VMR.xlsx --nt-fasta /path/to/viral.fna \
  --protein-fasta /path/to/viral.faa --plots supported
```
输出目录必须不存在。--plots 可取 supported（A/B）、all 或 none。

## 新组装数据
```bash
PYTHON_BIN=/path/to/python bash screen_virus_res.sh \
  --query assembly.fasta --reads1 R1.fq.gz --reads2 R2.fq.gz \
  --out new_analysis --vmr /path/to/VMR.xlsx \
  --nt-fasta /path/to/viral.fna --protein-fasta /path/to/viral.faa \
  --nt-db /path/to/blast/viral --viral-db /path/to/viral.dmnd \
  --competitive-db /path/to/competitive.dmnd \
  --competitive-fasta /path/to/competitive.fasta.gz \
  --taxonomy-map /path/to/taxonomy.tsv --pfam /path/to/Pfam-A.hmm
```
仅有组装时用 --no-mapping 替代两个 reads 参数。--preflight 仅检查发现阶段。环境变量 BLASTN_BIN、DIAMOND_BIN、HMMSEARCH_BIN、BWA_BIN、SAMTOOLS_BIN 可覆盖工具路径。数据库与分类映射不随源码发布，分类映射必须对应竞争库；不提供映射时可能大量 unclassified。

增强内核仍保留原主机的默认路径；跨机器必须显式配置，尚不是完全可移植发行版。VMR 解析仍沿用传统模块，对工作表名 VMR MSL41 有要求。

## 证据等级
A：至少两类 nt/protein/RdRP 证据、竞争最佳分类 viral，且 MAPQ20 覆盖广度>=0.8、平均深度>=5。
B：至少两类证据但不满足 A。
C：单类/待解决证据。
D：竞争最佳命中 cellular，优先标记冲突。
所有候选保留。分级属于启发式优先级，不是感染、物种或新病毒判定；同源性证据不完全独立。竞争证据只取最高分分类，尚未实现区域级分数竞争。VMR 的宿主信息仅描述参考病毒，不能证明当前植物感染。

## 输出
report/annotated_candidates.tsv：全部候选及原因、回贴和ICTV字段。
report/virus_evidence_summary.tsv：参考相关名称分组和等级计数，分组数不等于物种数。
report/virus_summary_legacy.tsv：兼容传统汇总。
report/nt_plots 和 protein_plots：HTML/SVG。
report/validation.json：计数、输入哈希和阈值。
run_manifest.json：实际运行状态；REPORT.md：解读边界。

## 验证范围
已用现有 LGTZ/HZS 完整候选表验证后处理；不重跑比对。新组装端到端执行和PDF导出尚未验证。validation_LGTZ 是绘图适配修复前的开发结果，以 validation_LGTZ_fixed 及后续 final 目录为准。
