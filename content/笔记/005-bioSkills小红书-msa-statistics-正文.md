<!--
META
用途: 005 msa-statistics 小红书帖子「正文文本框」文案，与出图源稿配套，不进站点、不参与 md2card。
标题建议: 蛋白IC别套DNA公式，差出2.8倍
/ META
-->

做多序列比对（MSA）之后，常要算每一列的信息量，看哪些位置保守、在进化上重要。bioSkills 的 msa-statistics 就是干这个的。

拿示例数据试了下，踩到一个点：

算蛋白序列每列的信息量（IC），不能直接套用 DNA 的均匀假设。氨基酸频率差很大（Leu 常见、Trp 稀有），uniform 把它们当一样——可稀有的残基保守住，信号本该更强，这层恰恰被 uniform 抹平了。

换成 Robinson 经验背景后，同一个保守位点：
Leu（常见残基）→ 3.44 bits
Trp（稀有残基）→ 6.23 bits
同样保守，稀有的信息量更高，差异这才显出来。

所以跑蛋白 MSA 的 IC 统计，背景用 Robinson，别用 uniform。

#生信 #生物信息学 #多序列比对 #MSA #bioSkills
