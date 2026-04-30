# 38 Targeted Enhancement Summary (Stable Mode)

## 本轮策略
- 采用分层精算：先 3 candidates × 2 endpoints（OS/PFI）。
- Cox 仅在通过稳定性门控后执行；否则跳过并记录原因。
- 不新增模型家族，不扩展 ge3 主线。

## 生存证据口径
- 主证据：KM + log-rank。
- Cox：辅助证据，仅在可收敛且不过分离时报告。
- DSS/DFI 状态：已在稳定条件下补算。

## 关键输出
- results/tables/multi_endpoint_clinical_validation_stable.csv
- results/tables/targeted_enhancement_survival_summary.csv
- results/figures/targeted_km_os_real.png
- results/figures/targeted_km_pfi_real.png
- results/figures/multi_endpoint_validation_comparison_stable.png

## 结论
- partial_ge4 未在稳定生存主证据（OS/PFI）上超过 old final baseline。
- 是否替换 old final baseline：仅当稳定主证据与其他维度形成一致优势时才考虑替换。