# 11 数据读取完整性校验

本文件基于 results/logs 中的结构化日志自动生成。

- 说明："头部预览"仅用于识别文件结构，不代表已完成全量读取。
- 说明："已完成全量扫描"表示该文件已通过流式方式完整遍历，维度统计可作为后续文档依据。
- 说明：若某文件尚未完成完整解析，当前维度仅为初步统计，待后续预处理阶段进一步确认。

## 1. 原始文件读取完整性总表

| key | file_size_mb | delimiter_guess | encoding_guess | total_lines | header_column_count | parsed_column_count | read_mode | full_scan_completed | notes |
|---|---:|---|---|---:|---:|---:|---|---|---|
| mutation | 34.242 | tab | utf-8-sig | 40544 | 440 | 440 | streamed_scan | True | streamed full scan completed |
| cnv | 22.997 | tab | utf-8-sig | 24777 | 442 | 442 | streamed_scan | True | streamed full scan completed |
| rna | 63.689 | tab | utf-8-sig | 20531 | 451 | 451 | streamed_scan | True | streamed full scan completed |
| methylation | 1158.537 | tab | utf-8-sig | 485578 | 399 | 399 | streamed_scan | True | large_file_stream_scan_recommended |
| methylation_probe_map | 17.361 | tab | utf-8-sig | 395986 | 6 | 6 | streamed_scan | True | streamed full scan completed |
| mirna | 5.577 | tab | utf-8-sig | 2179 | 429 | 429 | streamed_scan | True | streamed full scan completed |
| rppa | 1.05 | tab | utf-8-sig | 228 | 358 | 358 | streamed_scan | True | streamed full scan completed |
| clinical_matrix | 0.501 | tab | utf-8-sig | 581 | 108 | 108 | streamed_scan | True | streamed full scan completed |
| survival | 0.025 | tab | utf-8-sig | 512 | 11 | 11 | streamed_scan | True | streamed full scan completed |

## 2. 文件复制一致性校验

| key | size_match | hash_match | hash_algorithm | notes |
|---|---|---|---|---|
| mutation | True | True | sha256 | copy_verified |
| cnv | True | True | sha256 | copy_verified |
| rna | True | True | sha256 | copy_verified |
| methylation | True | True | sha256 | copy_verified |
| methylation_probe_map | True | True | sha256 | copy_verified |
| mirna | True | True | sha256 | copy_verified |
| rppa | True | True | sha256 | copy_verified |
| clinical_matrix | True | True | sha256 | copy_verified |
| survival | True | True | sha256 | copy_verified |

## 3. 风险点与后续建议

- 对超大矩阵（如 methylation）避免默认一次性全量 read_csv，优先使用 chunked_read 或 streamed_scan。
- 后续预处理脚本需记录输入维度、读取模式（full/chunked）、输出维度、过滤数量与原因。
- 若出现列数不一致，需在预处理阶段先做格式修复再进入建模。
- 可评估将超大中间结果转存为 parquet 以降低重复 I/O 成本。
