<!--
META
用途: 006 alignment-trimming 小红书帖子「正文文本框」文案，与出图源稿配套，不进站点、不参与 md2card。
标题建议: 比对修剪 ClipKIT 的默认行为
/ META
-->

多序列比对（MSA）建完后，常要去掉不可靠的列（gap 太多、保守度太低），让下游系统发育树更稳。bioSkills 的 alignment-trimming 就是干这个的。

拿仓库自带的示例 MSA，用 ClipKIT 默认推荐模式 kpic-smart-gap 试了下，记一个观察：

21 列只砍掉 1 列（第 13 列 gap 比例 0.75），保留率 95.2%——落在"安全修剪区"，对树几乎没影响。加 --log 会输出每列 kept/trim + gap 比例的审计日志；不加只拿到 trimmed.fasta，别人没法核对你到底砍了哪列。

所以比对修剪"轻修保信号、重修丢信号"，默认用 kpic-smart-gap，保留率掉到 70% 以下就换更温和的模式（或干脆不修）。

#生信 #生物信息学 #多序列比对 #MSA #ClipKIT #bioSkills
