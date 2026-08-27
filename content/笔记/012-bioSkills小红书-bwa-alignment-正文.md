# 012｜bioSkills bwa-alignment：DNA 比对前先注入 read group

<!--
META
用途: 012 bwa-alignment 小红书帖子「正文文本框」文案，与出图源稿配套，不进站点、不参与 md2card。
标题建议: bwa 比对别忘 read group
/META
-->

做 DNA 短读长比对，bwa-mem2 是变异调用流程的默认比对器。它最关键的不是「跑通」，而是比对前就把下游要用的契约写好——否则后续 GATK 直接报错或样本被静默混在一起。

实测三个发现：

① read group 是硬契约，不是装饰。比对时用 -R 注入 @RG（含样本名 SM、误差单元 ID、平台 PL、文库 LB），生成的 BAM 头里就有 @RG；不写，BAM 头里就完全没有 @RG。GATK 靠 SM 分组样本、靠 ID/LB 建模误差和去重，缺了要么报错要么行为错误，事后用 Picard 补是一次完整 BAM 重写。

② MAPQ 标度到 60，和 bowtie2 不一样。同一批数据实测 bwa-mem2 的 max MAPQ = 60，而 bowtie2 上限只有 42（end-to-end）/ 44（local）。所以「MAPQ ≥ 60 唯一映射」过滤在 bwa 上合理，照搬去压 bowtie2 的 BAM 会清空——阈值要跟着比对器走。

③ 去重顺序错了会静默标错。正确顺序是 collate（按名）→ fixmate -m（写 MC 标签）→ sort（坐标）→ markdup。fixmate 不带 -m 或拿坐标序 BAM 直接 markdup，重复标记会悄悄出错；SV 流程还要用 -Y 保留分裂读长全程序列，别用 -M（会把 SV 证据降级藏起来）。

结论：bwa 比对 = 为下游变异调用写对契约。比对即注入 @RG；人类数据用 GRCh38 + decoy；去重严格走 collate→fixmate -m→sort→markdup；SV 用 -Y 不用 -M；要可复现加 -K 100000000。
