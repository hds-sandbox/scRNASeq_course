# Introduction

scRNAseq data analysis of [Deng et
al. (2014)](https://www.science.org/doi/10.1126/science.1245316). This
study follows the early embryonic development of mouse cells from zygote
to late blastocists (and some adult cells).

![Mouse Early Development
<https://doi.org/10.1093/molehr/gaw002>](../../Notes/mouse_development.jpeg)

# Load libraries

    library(rhdf5)
    library(tidyverse)
    library(Seurat)
    library(SeuratWrappers)
    library(patchwork)
    library(biomaRt)
    library(slingshot)
    library(gprofiler2)

# Load dataset

Data is a matrix of cells in columns and genes in rows. We specify that
the row names (our genes) are in the first column.

    raw_data <- read.delim("../../Data/example_data/GSE45719_cts.txt", sep = "\t", row.names = 1)
    head(raw_data, n = 5)

    ##        GSM1112619 GSM1278045 GSM1112568 GSM1112558 GSM1112560 GSM1112735
    ## Hvcn1         427        363          2          0         37         36
    ## Gbp7            0          0          0          0          0          0
    ## Arrdc1        716          0       1414        920        446          0
    ## Ercc5           0       1918        283        190         29        264
    ## Mrpl15       3967        794       5338       5630       2941       4939
    ##        GSM1112655 GSM1278024 GSM1278023 GSM1112769 GSM1278019 GSM1278035
    ## Hvcn1        1266         42        148       1600         19          4
    ## Gbp7         1285          0          0        533          0          0
    ## Arrdc1       7679        431        411       1478        797         43
    ## Ercc5         166         33        116        779          0          0
    ## Mrpl15       1804       1083       1029         49       1856         81
    ##        GSM1112673 GSM1112606 GSM1112519 GSM1112720 GSM1112596 GSM1112557
    ## Hvcn1         168        966          0          1        801          4
    ## Gbp7            0       1301          0          0        102          0
    ## Arrdc1        173        305        374          0       3312        128
    ## Ercc5           0        463          3          0        109          0
    ## Mrpl15        326          0       6339       3799         36        469
    ##        GSM1112613 GSM1112535 GSM1112762 GSM1112711 GSM1112636 GSM1112573
    ## Hvcn1           1          0          1          5        648          0
    ## Gbp7            0          0          0          0          0          1
    ## Arrdc1        645       2374        699        102          0        680
    ## Ercc5          40          0       1904        102         53        257
    ## Mrpl15       2160       3868       2216       1244       4819       1790
    ##        GSM1112708 GSM1112659 GSM1112648 GSM1278037 GSM1112697 GSM1112513
    ## Hvcn1           2       1360          1          3        679          0
    ## Gbp7            0       1059          0         47        214          0
    ## Arrdc1          4      11805          0        255         97         42
    ## Ercc5           0        650          0          0         72          0
    ## Mrpl15       5210       3826        822        867        104       1988
    ##        GSM1112574 GSM1112641 GSM1112724 GSM1112651 GSM1112540 GSM1112645
    ## Hvcn1         166         18          0          2          0        104
    ## Gbp7            0          0          0          0          0          0
    ## Arrdc1        418          0          0          2       3677          0
    ## Ercc5         142          0       1328         28          0          0
    ## Mrpl15       4021       2771       5415       4083       2207       2723
    ##        GSM1112597 GSM1112672 GSM1112577 GSM1112588 GSM1278016 GSM1112500
    ## Hvcn1        1054          0        169          0          0          0
    ## Gbp7          201          0       2120         13          0          0
    ## Arrdc1       3531          0       1779          3        629        729
    ## Ercc5         179          0        214          7         27          0
    ## Mrpl15         21       1182       2743         64        921       2497
    ##        GSM1112550 GSM1112706 GSM1112490 GSM1112747 GSM1112666 GSM1112687
    ## Hvcn1         251          1          2        330          1          0
    ## Gbp7            0          0          0          0          0          0
    ## Arrdc1       3505          0       2246          0          0          0
    ## Ercc5           2          0          1          2          1          4
    ## Mrpl15        436       2010       2259       2550       3288         12
    ##        GSM1112682 GSM1112755 GSM1112661 GSM1112635 GSM1112494 GSM1112631
    ## Hvcn1           1          0       3337          0          0          0
    ## Gbp7            0          0        640          0          0          0
    ## Arrdc1        639         11        258        247       1682          2
    ## Ercc5         175          0        456          0        132          1
    ## Mrpl15       1008       2155        164       4268       5493       2583
    ##        GSM1112660 GSM1112594 GSM1112699 GSM1112571 GSM1112643 GSM1112549
    ## Hvcn1        1946          1        434          7        103        228
    ## Gbp7          361          0         45          1          0          8
    ## Arrdc1       1915          0        310        842          0       3580
    ## Ercc5           0         11          0        274        862          0
    ## Mrpl15       2523         33          0       3690       1719       3612
    ##        GSM1112679 GSM1112579 GSM1112756 GSM1278033 GSM1112683 GSM1112759
    ## Hvcn1           0        124          0          0          2          0
    ## Gbp7            0          0          0          0          0          0
    ## Arrdc1          0        843        219          0          0          3
    ## Ercc5           1         39        129          0          5          8
    ## Mrpl15        511        719        993        174       3864       4214
    ##        GSM1112497 GSM1112608 GSM1112758 GSM1112595 GSM1112703 GSM1112624
    ## Hvcn1           0       3525          1       1112        129        169
    ## Gbp7            0       1110          0         14         28          0
    ## Arrdc1       1776       1173          5       2323        331        223
    ## Ercc5         692        510        188        247         42          0
    ## Mrpl15       4769         53       1978         77          0       3000
    ##        GSM1112714 GSM1112640 GSM1112752 GSM1112495 GSM1112639 GSM1112518
    ## Hvcn1         386          0       1833          0          1          0
    ## Gbp7            0          0          0          0          0          0
    ## Arrdc1          0          0          0          0          0        932
    ## Ercc5          54          0          0          0        296          0
    ## Mrpl15       2320       2493       6517       5694       6254       6738
    ##        GSM1112617 GSM1112520 GSM1112561 GSM1112548 GSM1112665 GSM1278044
    ## Hvcn1           0         68       1234        110         47          0
    ## Gbp7            0          0          0          0          0        103
    ## Arrdc1          0        974       1537       3435          1          1
    ## Ercc5           0        234          0         12         59        314
    ## Mrpl15       5163       3922       2994       1714       1929        270
    ##        GSM1278040 GSM1278021 GSM1112610 GSM1112504 GSM1112709 GSM1112615
    ## Hvcn1           0         80       2352        289          0          1
    ## Gbp7            0          0        594          0          0          0
    ## Arrdc1          0        747        920        362          5          2
    ## Ercc5         207         18       1304         84        338          0
    ## Mrpl15        733        797          0       2821       3088       1669
    ##        GSM1112732 GSM1112725 GSM1112657 GSM1112632 GSM1112562 GSM1112633
    ## Hvcn1         934          0        398          0        303          1
    ## Gbp7            0          0        627          0          0          0
    ## Arrdc1          0          0        938          0       2112          0
    ## Ercc5           0          1          9        360          0        529
    ## Mrpl15       4970       1550       1697       1665       2546       1008
    ##        GSM1112700 GSM1112715 GSM1112618 GSM1112765 GSM1112623 GSM1278028
    ## Hvcn1         200         20          1          0         26         90
    ## Gbp7           73          0          0          0          0          0
    ## Arrdc1        357          0        402          3          8       2010
    ## Ercc5          39       1151         80          5         11          1
    ## Mrpl15         38        121       1434       2866       8078       1043
    ##        GSM1112537 GSM1112719 GSM1112547 GSM1112605 GSM1278020 GSM1112689
    ## Hvcn1           1          0        253       1294         37         20
    ## Gbp7            0          0          0        442          0          0
    ## Arrdc1        815          0       3356        511        263          0
    ## Ercc5           0          6          0        264          0          0
    ## Mrpl15       2634          0       1853          0        506       2677
    ##        GSM1112738 GSM1112701 GSM1278022 GSM1112707 GSM1112750 GSM1112509
    ## Hvcn1          96        314          1        209          0          0
    ## Gbp7            1          6          0          0          0          0
    ## Arrdc1        406        374        809        737          6       1939
    ## Ercc5          24         86          0          5        105          0
    ## Mrpl15       1640          0       1048       1105       3036       8395
    ##        GSM1112565 GSM1112531 GSM1112628 GSM1112525 GSM1112609 GSM1112602
    ## Hvcn1          11          2          2          0       2369        565
    ## Gbp7            0          0          0          0        518         61
    ## Arrdc1       1795         28          0       1070        733       1081
    ## Ercc5           0          0          1        109        501          0
    ## Mrpl15       2737       8015       5677       4656         29          0
    ##        GSM1278009 GSM1278042 GSM1112511 GSM1112702 GSM1112582 GSM1112589
    ## Hvcn1          72          0          1        251          0          0
    ## Gbp7            0        530          0         79        356          9
    ## Arrdc1          0          2        922        272          2          1
    ## Ercc5           0        330          0          2        303          4
    ## Mrpl15       5907       2083       6017          0       1737         67
    ##        GSM1112650 GSM1112730 GSM1112534 GSM1112522 GSM1112629 GSM1112642
    ## Hvcn1           1         51          0          2          1        693
    ## Gbp7            0         34          0          0          0          0
    ## Arrdc1          0        193          0          0          0          0
    ## Ercc5        5436        140          0          0        228          0
    ## Mrpl15       5310       2696       6200       4822       1991       2265
    ##        GSM1112492 GSM1112695 GSM1112521 GSM1112693 GSM1112634 GSM1112542
    ## Hvcn1           0        407          0          0          0         19
    ## Gbp7            0        126          4          0          0          0
    ## Arrdc1       1410         77        334         13        687       2651
    ## Ercc5           0         55          0         39          0        163
    ## Mrpl15       3969        107       4100       2268       2255       1266
    ##        GSM1112586 GSM1112646 GSM1112670 GSM1112587 GSM1278026 GSM1112498
    ## Hvcn1           9       1205         13         25        266        125
    ## Gbp7          801          0          0        441          0          0
    ## Arrdc1         15          0         26         17        920       1547
    ## Ercc5          76         67         25        187          0        717
    ## Mrpl15       3704       2872        298       2593       2469       3845
    ##        GSM1112675 GSM1112593 GSM1112512 GSM1278036 GSM1112545 GSM1112499
    ## Hvcn1           1          2         11          1        165          0
    ## Gbp7            0          0          1         10          0          0
    ## Arrdc1        433          0        476          3       1240        618
    ## Ercc5         385         29          0       1203          6        918
    ## Mrpl15       1877         19       4660        957       3102       1681
    ##        GSM1112739 GSM1112505 GSM1112754 GSM1278018 GSM1112612 GSM1112553
    ## Hvcn1         662          0          1        170          1          0
    ## Gbp7            3          0          0          0          0        553
    ## Arrdc1          2          5          0        115          0       3749
    ## Ercc5         325          0        174          0          0          0
    ## Mrpl15       4276      10683       3826       1239          1       2069
    ##        GSM1112572 GSM1112567 GSM1278010 GSM1112538 GSM1112740 GSM1112751
    ## Hvcn1           0         33          0        226         90          0
    ## Gbp7            0         18          0          0          0          0
    ## Arrdc1       2051       2074        180        213          5       1777
    ## Ercc5           0          0          0          0        161          0
    ## Mrpl15       5692       6095       3098       2309       2390       2585
    ##        GSM1112721 GSM1112625 GSM1278013 GSM1112546 GSM1112749 GSM1112591
    ## Hvcn1         297          8         21          0          0         29
    ## Gbp7            0          1          0        452          0          0
    ## Arrdc1          0          0          0       5050          0          2
    ## Ercc5           0          0          1        106          0          2
    ## Mrpl15       1253        930        873       3317       2031        148
    ##        GSM1278031 GSM1112622 GSM1112698 GSM1112662 GSM1112528 GSM1112680
    ## Hvcn1          88        195        313       2907       1302         97
    ## Gbp7            0          0        274          0          0          0
    ## Arrdc1       1969         24        105       7222        765          0
    ## Ercc5           0          0         76        258          0          0
    ## Mrpl15       2858       3806        135       1750       5053       3299
    ##        GSM1112722 GSM1112704 GSM1112599 GSM1112658 GSM1112530 GSM1278034
    ## Hvcn1           0        140       1018        775          0        116
    ## Gbp7            0         60        208        441          0          0
    ## Arrdc1          0        127       2633       2323        232        120
    ## Ercc5           0         13        259          0          0         95
    ## Mrpl15       2162         13         32          0       8495        317
    ##        GSM1112569 GSM1112529 GSM1112576 GSM1112696 GSM1112523 GSM1112767
    ## Hvcn1           0        159         17        467        368       2995
    ## Gbp7            0          0          0        206          0        473
    ## Arrdc1       3848        293       1682         92       1106       2481
    ## Ercc5         559          0         94          0        227       1549
    ## Mrpl15       8370       1210       2832         85       9447          0
    ##        GSM1112713 GSM1278014 GSM1112668 GSM1112502 GSM1112718 GSM1278030
    ## Hvcn1         368          0          0          0        140          0
    ## Gbp7            0          0          0          0          0          0
    ## Arrdc1          0          0        378        656         24        857
    ## Ercc5           2          0         57        216        323          0
    ## Mrpl15       1709       1129        513       2524       3681       2561
    ##        GSM1112681 GSM1278029 GSM1112627 GSM1112647 GSM1112524 GSM1112585
    ## Hvcn1         155          0          0         66        130          0
    ## Gbp7            0          0          0          0          0          0
    ## Arrdc1        515       1017          0         39       1177         31
    ## Ercc5          34          2          0         21          0          0
    ## Mrpl15       2263        622       2096        535       7365        542
    ##        GSM1112746 GSM1112491 GSM1112551 GSM1112532 GSM1112496 GSM1112563
    ## Hvcn1           0        147         13         60          0        108
    ## Gbp7            0          0         42          0          0          0
    ## Arrdc1          0        746       2355        630       1284       3471
    ## Ercc5           1          1        317          0          0        251
    ## Mrpl15          8       1614        985       1520       3783       5833
    ##        GSM1112630 GSM1112638 GSM1278038 GSM1278015 GSM1112686 GSM1112541
    ## Hvcn1           0          0          0        376         59        121
    ## Gbp7            0          0        880          0          0          0
    ## Arrdc1        140          0        910        415          6       1679
    ## Ercc5           0          0         44        208       1027        431
    ## Mrpl15       3092       3176        892       1193       1536       2030
    ##        GSM1112616 GSM1112694 GSM1278012 GSM1112539 GSM1112652 GSM1112644
    ## Hvcn1         755        227        161          0          1        877
    ## Gbp7            0        285          0          0          0          0
    ## Arrdc1          0        191        430          0         45          0
    ## Ercc5         597         32          1          0        236        436
    ## Mrpl15       4629        152       3110       1276       1889       2347
    ##        GSM1112712 GSM1112678 GSM1112729 GSM1112649 GSM1112705 GSM1112763
    ## Hvcn1           2          1          1          0        379          0
    ## Gbp7            0          0          0          0         50          0
    ## Arrdc1          0          0          0          0        907         12
    ## Ercc5         126          0        117          2        177          0
    ## Mrpl15       1504       1746       2253         34         14       2286
    ##        GSM1112743 GSM1112570 GSM1278032 GSM1112717 GSM1112690 GSM1112663
    ## Hvcn1           0          0        106          0        264       1037
    ## Gbp7            0          0          0          0          1        336
    ## Arrdc1        198       1928          4          1          1       6839
    ## Ercc5         101        180          1        198        791       1040
    ## Mrpl15        827       6005         30       2816       1323       5044
    ##        GSM1112671 GSM1112757 GSM1112566 GSM1278043 GSM1112764 GSM1112742
    ## Hvcn1         144          0          0          1         39        131
    ## Gbp7            0          0          0        474          0          0
    ## Arrdc1         44        545         70        354        454         34
    ## Ercc5           0          0          0          1        610          3
    ## Mrpl15       2010        476       2289        902       1605       2043
    ##        GSM1112543 GSM1112626 GSM1112637 GSM1112736 GSM1112748 GSM1278011
    ## Hvcn1           2         80          2        381          0          0
    ## Gbp7           17          0          0          0          0          0
    ## Arrdc1       1223          2        634          0          1        485
    ## Ercc5          69        387        115          0          0          0
    ## Mrpl15       2311       2403       2208        678         39       3846
    ##        GSM1112710 GSM1112723 GSM1112691 GSM1112737 GSM1112515 GSM1112669
    ## Hvcn1           4          0          6          0          0         12
    ## Gbp7            0          0          0          0          0          0
    ## Arrdc1          0          0          7        632       1565         67
    ## Ercc5           3         14          0          0          0          0
    ## Mrpl15       3355       1749       1362       6766       4536        867
    ##        GSM1112674 GSM1112685 GSM1112575 GSM1112734 GSM1112503 GSM1112676
    ## Hvcn1         402          9        234          0        214         85
    ## Gbp7            0          0          0          0          0          0
    ## Arrdc1          0          1          6        157        663          0
    ## Ercc5         360        244          0         90          0        606
    ## Mrpl15       1708       3834       1218       3616       2371       3498
    ##        GSM1112684 GSM1112517 GSM1112592 GSM1112556 GSM1112544 GSM1112555
    ## Hvcn1           1          0          0          8         36        386
    ## Gbp7            0          0          2          1          0          0
    ## Arrdc1          2         17        110        144       1265        445
    ## Ercc5           0        127          0          0          6        150
    ## Mrpl15       1746       4899         61        612       1006       2943
    ##        GSM1112768 GSM1112590 GSM1112583 GSM1278039 GSM1112598 GSM1112664
    ## Hvcn1        2290          2          0        194        477          1
    ## Gbp7          719          0         34          0         81          0
    ## Arrdc1       2269          2          0        662       1954          1
    ## Ercc5        1601          1         14         69        149          2
    ## Mrpl15         42          9        217        755        168       2727
    ##        GSM1112533 GSM1112692 GSM1112584 GSM1112580 GSM1112552 GSM1112760
    ## Hvcn1         270        292          0        745        391          0
    ## Gbp7            1          0          1          0          1          0
    ## Arrdc1       1123        770          0       1304       2119        317
    ## Ercc5           0          0          0         44          2          0
    ## Mrpl15       4842       2828        273       3120       1585       2884
    ##        GSM1112656 GSM1112600 GSM1112766 GSM1112677 GSM1112578 GSM1112744
    ## Hvcn1        1599        726       3204         13         87          0
    ## Gbp7            1         84       1101          0          1          0
    ## Arrdc1       8422       1749       3574          0        321         55
    ## Ercc5           7        136       1391          0          4          1
    ## Mrpl15       1501        114          0         93       1389       1169
    ##        GSM1112601 GSM1112716 GSM1112604 GSM1278017 GSM1112654 GSM1112603
    ## Hvcn1         538          0       1048        134       1304       1272
    ## Gbp7            0          0        325          0        392        261
    ## Arrdc1       1669          0        176        917       1076        235
    ## Ercc5           0         21         53          1        529        486
    ## Mrpl15         47       1459          0       1414       1548          3
    ##        GSM1112653 GSM1112607 GSM1112611 GSM1112516 GSM1112506 GSM1278027
    ## Hvcn1         296       1248        157          0          0         90
    ## Gbp7            0       1028          0          0          0          0
    ## Arrdc1        617        673          4          0        399        725
    ## Ercc5          14        114        762        719          0          0
    ## Mrpl15       1035        128       1928       7617       1188       1519
    ##        GSM1112621 GSM1112614 GSM1112514 GSM1112726 GSM1112745 GSM1112501
    ## Hvcn1         792          0          0          0          0          0
    ## Gbp7            0          0          0          0          1          0
    ## Arrdc1          1          2          1          0          3       1151
    ## Ercc5           1          2          0          0          0        350
    ## Mrpl15       5086        129        518       3181        288       1376
    ##        GSM1112667 GSM1112508 GSM1112733 GSM1278025 GSM1112564 GSM1112493
    ## Hvcn1           0          0        279          0          0        103
    ## Gbp7            0          0          0          0          0          0
    ## Arrdc1          1        315          0        254        359        727
    ## Ercc5           0          0         12          0          2          0
    ## Mrpl15        219       4234       2734       1168       5881       3657
    ##        GSM1112753 GSM1112554 GSM1112728 GSM1112727 GSM1112741 GSM1112527
    ## Hvcn1           7          6         23        151          1          0
    ## Gbp7            1          0          0          0          0          0
    ## Arrdc1          0        827          0          0          0          0
    ## Ercc5           1       1298         92       1421          0        295
    ## Mrpl15       2106       4400       1772       2949       1504       2979
    ##        GSM1112507 GSM1112688 GSM1112559 GSM1112761 GSM1112731 GSM1112620
    ## Hvcn1        1541          0         37         17        560        400
    ## Gbp7           50          0          0          0          0          0
    ## Arrdc1          0          0        303       1412        230        110
    ## Ercc5           0          0        189          0        328        148
    ## Mrpl15          0       4815       4404       2541       2147       2780
    ##        GSM1112526 GSM1112510 GSM1278041 GSM1112536 GSM1112581
    ## Hvcn1           0        863          0         55         65
    ## Gbp7            0          0          0          0          0
    ## Arrdc1          0        984          0        380       1031
    ## Ercc5           6        289       1026        332          0
    ## Mrpl15       5308       9164       1952       4788       2260

Here is the metadata. The column names of `raw_data` are the same as the
row names of metadata.

    metadata <- read.delim("../../Data/example_data/GSE45719_metadata.tsv", header = T, sep = "\t", quote = "", stringsAsFactors = F, row.names = 1)
    head(metadata, n = 5)

    ##                             Strain                  Cross        Cell_type
    ## GSM1112619 CAST/EiJ(m)xC57BL/6J(f) first generation cross Early blastocyst
    ## GSM1278045 CAST/EiJ(m)xC57BL/6J(f) first generation cross       Fibroblast
    ## GSM1112568 CAST/EiJ(m)xC57BL/6J(f) first generation cross           8-cell
    ## GSM1112558 CAST/EiJ(m)xC57BL/6J(f) first generation cross           8-cell
    ## GSM1112560 CAST/EiJ(m)xC57BL/6J(f) first generation cross           8-cell

For the purposes of this workshop, the metadata of this project has been
simplified. The most interesting variable is **Cell\_type**, which shows
the developmental stage of our cells. Cell type is not ordered by its
differentiation stage, but we can do that now:

    metadata$Cell_type <- factor(metadata$Cell_type, 
                                 levels = c("MII Oocyte","Zygote",
                                            "Early 2-cell","Mid 2-cell","Late 2-cell",
                                            "4-cell","8-cell","16-cell",
                                            "Early blastocyst","Mid blastocyst","Late blastocyst", 
                                            "Fibroblast","Adult"))

There are some cell types that we are not really interested in, and will
be removed in the next step below.

## Seurat object

Now that we have our data loaded, we can create a **Seurat Object**,
from which we will perform our analysis. The function
`CreateSeuratObject()` takes as arguments our count matrix, metadata and
several options for filtering the data:

-   `min.cells` will keep features (genes) that are detected in at least
    this many cells.
-   `min.features` will keep cells with at least this many features
    detected.

This will prevent from keeping cells and genes with an immense majority
of 0’s as values. In addition, we are removing some cell types we do not
want to analyse.

    raw_ann <- CreateSeuratObject(counts = raw_data, meta.data = metadata, min.cells = 3, min.features = 200,)
    raw_ann <- subset(raw_ann, cells = colnames(raw_ann)[!raw_ann$Cell_type %in% c("MII Oocyte","Fibroblast","Adult")])

# Quality Control

The Seurat object initialization step above only considered cells that
expressed at least 300 genes and genes detected in at least 3 cells.
Here is how many cells and genes we start with:

    print(paste0("Before filtering: ", dim(raw_ann)[2], " cells ",  dim(raw_ann)[1], " genes"))

    ## [1] "Before filtering: 286 cells 20566 genes"

Additionally, we would like to exclude cells that are damaged. A common
metric to judge this (although by no means the only one) is the relative
expression of mitochondrially derived genes. When the cells apoptose due
to stress, their mitochondria becomes leaky and there is widespread RNA
degradation. Thus a relative enrichment of mitochondrially derived genes
can be a tell-tale sign of cell stress. Here, we compute the proportion
of transcripts that are of mitochondrial origin for every cell
(*percent.mito*), and visualize its distribution as a violin plot. We
also use the `GenePlot()` function to observe how percent.mito
correlates with other metrics.

    raw_ann[['percent.mito']] <- PercentageFeatureSet(raw_ann, pattern = "^mt-")

Now, sometimes, mitochondrial genes are a bit tricky to find, specially
if your genes are not gene names, but gene IDs. You might have already a
collection of gene names you want to use, but it is not always the case.
Therefore, it might very useful to have some extra annotation on our
genes which will help to select mitochondrial genes (or ERCC genes and
ribosomal genes). Since this dataset was aligned using the mouse genome
version *mm9*, we will use the mm9 annotation from
[Biomart](http://www.biomart.org/). The package **biomaRt** will do the
job!

    ensembl67=useMart(host='http://feb2014.archive.ensembl.org/', # select latest version of mm9 genome annotation
                      biomart='ENSEMBL_MART_ENSEMBL', dataset = "mmusculus_gene_ensembl")

    # we get a list of annotations we would like to fetch
    mm9.gene.annotations <- biomaRt::getBM(mart = ensembl67, attributes=c("ensembl_gene_id", "external_gene_id", "description", "chromosome_name"))
    head(mm9.gene.annotations)

    ##      ensembl_gene_id external_gene_id
    ## 1 ENSMUSG00000095309         Vmn1r125
    ## 2 ENSMUSG00000000126            Wnt9a
    ## 3 ENSMUSG00000086196          Gm13571
    ## 4 ENSMUSG00000054418    2900041M22Rik
    ## 5 ENSMUSG00000095268           Gm2913
    ## 6 ENSMUSG00000082399          Gm14036
    ##                                                                  description
    ## 1             vomeronasal 1 receptor 125 [Source:MGI Symbol;Acc:MGI:3704112]
    ## 2 wingless-type MMTV integration site 9A [Source:MGI Symbol;Acc:MGI:2446084]
    ## 3                   predicted gene 13571 [Source:MGI Symbol;Acc:MGI:3652223]
    ## 4             RIKEN cDNA 2900041M22 gene [Source:MGI Symbol;Acc:MGI:1925653]
    ## 5                    predicted gene 2913 [Source:MGI Symbol;Acc:MGI:3781091]
    ## 6                   predicted gene 14036 [Source:MGI Symbol;Acc:MGI:3649581]
    ##   chromosome_name
    ## 1               7
    ## 2              11
    ## 3               2
    ## 4              11
    ## 5               X
    ## 6               2

We select that mitochondrial genes are those with the
**chromosome\_name** annotation equal to “MT”.

    mt_genes <- mm9.gene.annotations %>% filter(chromosome_name == "MT") %>% pull(external_gene_id)
    mt_genes <- mt_genes %in% rownames(raw_ann)

    raw_ann[['percent.mito']] <- PercentageFeatureSet(raw_ann, features = mt_genes)

If your data contains spike-ins (ERCC genes), you can also compute the
proportion of transcripts belonging to them. Furthermore, it might be
useful to calculate the percentage of reads aligned to ribosomal genes,
since it has been shown that they can skew the data due to their high
variability.

    raw_ann[['percent.ercc']] <- PercentageFeatureSet(raw_ann, pattern = "^ERCC-")
    raw_ann[['percent.ribo']] <- PercentageFeatureSet(raw_ann, pattern = "^Rp[ls]")

Nonetheless, in this experiment there aren’t ERCC genes or mitochondrial
genes:

    sum(raw_ann$percent.ercc)

    ## [1] 0

    sum(raw_ann$percent.mito)

    ## [1] 0

## Visualizations

It is extremely useful to visualize QC measurements calculated so far.
Violin plots (fancy boxplots) are a great way to check the distribution
of values of all our QC measurements. When initializing the Seurat
Object, Seurat calculates also the number of genes detected and the
total library size per cell.

    VlnPlot(raw_ann, 
            features = c("nFeature_RNA", "nCount_RNA", "percent.ribo"),
            ncol = 4, group.by = "Cell_type")

<img src="scRNAseq_Seurat_files/figure-markdown_strict/QC_before_filter-1.png" style="display: block; margin: auto;" />

It seems that the late 2-cells have a very high total number of reads.
We should not filter them out, there are very few! In addition, our
percentage of ribosomal counts is also quite low (maximum is ~5%).

We can check the relationship between library size and number of genes
detected:

    FeatureScatter(raw_ann, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Cell_type")

<img src="scRNAseq_Seurat_files/figure-markdown_strict/QC_before_filtering2-1.png" style="display: block; margin: auto;" />
In other experiment with more cells, we might consider removing cells
with an unexpectedly high number of detected genes and reads. These
cells might be doublets, that is, two cells that were sequenced
together!

## Filtering

Some automatic filtering can be made using quantiles! Using the extreme
upper and lower quantiles (0.99 and 0.01) we can make sure that outliers
are removed. In this case, we will remove cells with very low library
size and genes detected (lower than quantile 0.01).

    feature_min <- quantile(raw_ann$nFeature_RNA, probs = 0.01)
    count_min <- quantile(raw_ann$nCount_RNA, probs = 0.01)

We can subset our dataset using the function `subset()`.

    adata <- subset(raw_ann, subset = 
                        nFeature_RNA > feature_min  & 
                        nCount_RNA > count_min)

    rm(raw_ann) # we remove the initial unfiltered dataset to reduce computational resources, this is not necessary!

Finally, this is how many cells and genes we have after filtering:

    print(paste0("After filtering: ", dim(adata)[2], " cells ",  dim(adata)[1], " genes"))

    ## [1] "After filtering: 280 cells 20566 genes"

And this is how our filtered data looks like:

    VlnPlot(adata, 
            features = c("nFeature_RNA", "nCount_RNA", "percent.ribo"),
            ncol = 4, group.by = "Cell_type")

<img src="scRNAseq_Seurat_files/figure-markdown_strict/QC_after_filtering-1.png" style="display: block; margin: auto;" />

# Exploratory analysis

## Normalization

Now that the data is filtered, we can proceed to normalize our count
matrix. Seurat normalizes the gene expression measurements for each cell
by the total expression, multiplies this by a scale factor (10,000 by
default), and log-transforms the result. There have been many methods to
normalize the data, but this is the simplest and the most intuitive. The
division by total expression is done to change all expression counts to
a relative measure, since experience has suggested that technical
factors (e.g. capture rate, efficiency of reverse transcription) are
largely responsible for the variation in the number of molecules per
cell, although genuine biological factors (e.g. cell cycle stage, cell
size) also play a smaller, but non-negligible role. The
log-transformation is a commonly used transformation that has many
desirable properties, such as variance stabilization (can you think of
others?).

    adata <- NormalizeData(adata)
    adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = 2000)

## Identification of Variable Genes

Identify most variable genes and label top 5 most highly variable.

    top10 <- head(VariableFeatures(adata), 5)
    LabelPoints(plot = VariableFeaturePlot(adata), points = top10, repel = T)

<img src="scRNAseq_Seurat_files/figure-markdown_strict/top_variable_genes-1.png" style="display: block; margin: auto;" />

## Scale

Gene expression scaling is necessary for proper clustering of our cells.
Since genes may be expressed in very different orders of magnitude,
extreme expression levels may drive the separation between cells and
bias the results.

    adata <- ScaleData(adata)

## PCA

We can perform a Principal Component analysis using the following
function. PCA can be used for visualization of our cells, as well as
clustering and other dimensionality reduction methods such as t-SNE or
UMAP.

    adata <- RunPCA(adata)

We can also calculate and visualize the importance of each gene for each
Principal Component.

    pcalod_1 <- VizDimLoadings(object = adata, dims = 1) + theme(axis.text.y = element_text(size = 8)) 
    pcalod_2 <- VizDimLoadings(object = adata, dims = 2) + theme(axis.text.y = element_text(size = 8))

    CombinePlots(plots = list(pcalod_1, pcalod_2), ncol = 2)

<img src="scRNAseq_Seurat_files/figure-markdown_strict/PCA_loadings-1.png" style="display: block; margin: auto;" />

Visualization only show negative loadings cause there are many more than
positives. We can see a balanced plot using *balance = TRUE*

    pcalod_1 <- VizDimLoadings(object = adata, dims = 1, balanced = TRUE) + theme(axis.text.y = element_text(size = 8)) 
    pcalod_2 <- VizDimLoadings(object = adata, dims = 2, balanced = TRUE) + theme(axis.text.y = element_text(size = 8))

    CombinePlots(plots = list(pcalod_1, pcalod_2), ncol = 2)

<img src="scRNAseq_Seurat_files/figure-markdown_strict/PCA_loadings_balanced-1.png" style="display: block; margin: auto;" />

## Cell clustering

Now that we have calculated our components, we can proceed to select the
number of PCs necessary to perform clustering. There are two methods
that we can use to determine the proper number of dimensions:

### PC selection

**Jackstraw method**

The Jackstraw method randomly permutes a subset of data, and calculates
projected PCA scores for these ‘random’ genes. Then compares the PCA
scores for the ‘random’ genes with the observed PCA scores to determine
statistical significance. End result is a p-value for each gene’s
association with each principal component.

    adata <- JackStraw(adata, num.replicate = 100)
    adata <- ScoreJackStraw(adata, dims = 1:20)
    JackStrawPlot(adata, dims = 1:20)

<img src="scRNAseq_Seurat_files/figure-markdown_strict/jackstraw-1.png" style="display: block; margin: auto;" />

As you can see in the plot above, there is a change in the orders of
magnitude of the PCs’ p-values. Suggesting that we may cut off around
PC8.

**Elbow method**

The elbow method allows us to explore the explained variation of each of
the Principal Components. The plot usually looks like an “elbow”, where
adding more PCs does not really contribute to the amount of explained
variation. We can see again that we reach a plateau around PC6 or PC8.

    ElbowPlot(adata)

<img src="scRNAseq_Seurat_files/figure-markdown_strict/elbow_plot-1.png" style="display: block; margin: auto;" />

### Clustering

We could include 6 PCs, but it does not hurt to use more; there are very
few cells in this experiment and we would like to include as much
information as possible. Use the `resolution` argument of the
`FindClusters()` function to fine-tune the number of clusters to find.
The larger the number the less clusters it will find.

    adata <- FindNeighbors(adata, dims = 1:20)
    adata <- FindClusters(adata, resolution = 0.8)

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 280
    ## Number of edges: 6771
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.7890
    ## Number of communities: 7
    ## Elapsed time: 0 seconds

# Visualization

    adata <- RunTSNE(adata)
    adata <- RunUMAP(adata, dims = 1:20)

## PCA

    DimPlot(adata, reduction = "pca", group.by = c('Cell_type')) + 
      DimPlot(adata, reduction = "pca", group.by = c('seurat_clusters'))

<img src="scRNAseq_Seurat_files/figure-markdown_strict/PCA_plot-1.png" style="display: block; margin: auto;" />

## UMAP

    p1 <- DimPlot(adata, reduction = "umap", group.by = 'Cell_type')
    p2 <- DimPlot(adata, reduction = "umap")
    p1 + p2

<img src="scRNAseq_Seurat_files/figure-markdown_strict/UMAP_plot-1.png" style="display: block; margin: auto;" />

## TSNE

    p3 <- DimPlot(adata, reduction = "tsne", group.by = 'Cell_type')
    p4 <- DimPlot(adata, reduction = "tsne")
    p3 + p4

<img src="scRNAseq_Seurat_files/figure-markdown_strict/TSNE_plot-1.png" style="display: block; margin: auto;" />

# Differential Expression Analysis

Differential testing works a bit different than in bulk RNA-Seq
analysis. Usually, your experiment will include dozens or hundreds of
cells per cluster/condition. We can make use of Wilcox tests and
multiple testing correction to identify statistically significant genes
using the `FindAllMarkers()` function. This function will gather all
“markers” for each of your *Ident* variable values (your conditions or
identified clusters).

    markers <- FindAllMarkers(adata,logfc.threshold = 0.5, only.pos = T)

    markers %>%
        group_by(cluster) %>%
        top_n(n = 10, wt = avg_log2FC) -> top10

    head(top10)

    ## # A tibble: 6 × 7
    ## # Groups:   cluster [1]
    ##      p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene 
    ##      <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>
    ## 1 8.95e-29       2.76     1 0.583  1.84e-24 0       Upp1 
    ## 2 1.44e-28       2.24     1 0.413  2.96e-24 0       Spp1 
    ## 3 2.30e-28       1.94     1 0.955  4.74e-24 0       Uhrf1
    ## 4 6.59e-28       2.42     1 0.798  1.36e-23 0       Tdgf1
    ## 5 1.65e-26       1.93     1 0.628  3.39e-22 0       Glrx 
    ## 6 4.63e-26       1.89     1 0.928  9.52e-22 0       Fabp5

On the other hand, if you want to make a specific comparison, you can
use the `FindMarkers()` function:

    cluster2.markers <- FindMarkers(adata, ident.1 = 2, min.pct = 0.5, only.pos = T) # Only cluster 2 markers
    cluster2.markers %>% arrange(desc(avg_log2FC)) %>% head(10)

    ##                p_val avg_log2FC pct.1 pct.2    p_val_adj
    ## Gm11517 1.107520e-16   1.755328     1 0.975 2.277726e-12
    ## Impdh2  2.399680e-21   1.683743     1 0.992 4.935181e-17
    ## Larp4   2.697206e-19   1.664780     1 0.987 5.547075e-15
    ## Odc1    3.235495e-22   1.621770     1 1.000 6.654119e-18
    ## Obox6   3.386860e-15   1.545603     1 0.574 6.965416e-11
    ## Ppt2    4.479082e-19   1.534751     1 0.772 9.211680e-15
    ## Gldc    1.116913e-20   1.522874     1 1.000 2.297044e-16
    ## Trim43b 9.055909e-18   1.478628     1 0.983 1.862438e-13
    ## Trim43c 5.894929e-18   1.474335     1 0.966 1.212351e-13
    ## Sp110   7.422250e-17   1.472109     1 0.945 1.526460e-12

Even compare specific clusters against others by selecting `ident.1` and
`ident.2`:

    cluster5.markers <- FindMarkers(adata, ident.1 = 5, ident.2 = c(0, 3), # Differences between cluster 5 and clusters 0 and 3
                                    min.pct = 0.5, only.pos = T) 

## Visualizations

We can visualize all markers in a heatmap just like this! To go easy on
the visualization, we only use the top 10 markers.

    DoHeatmap(adata, features = top10$gene) + NoLegend()

<img src="scRNAseq_Seurat_files/figure-markdown_strict/DEA_heatmap-1.png" style="display: block; margin: auto;" />

We can also plot the expression on our PCA, tSNE or UMAP. We only need
to pass a vector of genes to the `features` argument of `FeaturePlot()`

    FeaturePlot(adata, features = top10$gene[1:4], reduction = "pca")

<img src="scRNAseq_Seurat_files/figure-markdown_strict/DEA_feature-1.png" style="display: block; margin: auto;" />

Or as a violin plot, a ridge plot or a dot plot.

    VlnPlot(adata, features = top10$gene[1]) + RidgePlot(adata, features = top10$gene[1])

<img src="scRNAseq_Seurat_files/figure-markdown_strict/DEA_violin-1.png" style="display: block; margin: auto;" />

    DotPlot(adata, features = top10$gene[1:10])

<img src="scRNAseq_Seurat_files/figure-markdown_strict/DEA_dotplot-1.png" style="display: block; margin: auto;" />

# Functional analysis with gprofiler2

`gost()` function allows us to do functional profiling of gene lists,
such as our differentially expressed genes. The function performs
statistical enrichment analysis to find over-representation of terms
from Gene Ontology, biological pathways like KEGG and Reactome, human
disease annotations, etc. This is done by using hypergeometric tests
that are corrected for multiple testing.

## Single query

A standard input of the `gost()` function is a (named) list of gene
identifiers. The list can consist of mixed types of identifiers
(proteins, transcripts, microarray IDs, etc), SNP IDs, chromosomal
intervals or functional term IDs.

The result is a named list where *result* is a data.frame with the
enrichment analysis results and *meta* containing a named list with all
the metadata for the query.

    gostres <- gost(query = rownames(cluster2.markers), 
                    organism = "mmusculus", ordered_query = FALSE, 
                    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = FALSE, 
                    user_threshold = 0.05, correction_method = "g_SCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = NULL, as_short_link = FALSE)

    head(gostres$result)

    ##     query significant      p_value term_size query_size intersection_size
    ## 1 query_1        TRUE 5.266520e-15        46        100                26
    ## 2 query_1       FALSE 1.268504e-01         7        100                 4
    ## 3 query_1       FALSE 1.268504e-01         7        100                 4
    ## 4 query_1       FALSE 5.461155e-01         2        100                 2
    ## 5 query_1       FALSE 5.461155e-01         2        100                 2
    ## 6 query_1       FALSE 5.461155e-01         2        100                 2
    ##   precision    recall    term_id source                            term_name
    ## 1      0.26 0.5652174 CORUM:3047  CORUM Parvulin-associated pre-rRNP complex
    ## 2      0.04 0.5714286  CORUM:537  CORUM                     Mediator complex
    ## 3      0.04 0.5714286 CORUM:3029  CORUM                       Drosha complex
    ## 4      0.02 1.0000000 CORUM:2545  CORUM                   Grb2-mSos1 complex
    ## 5      0.02 1.0000000 CORUM:2186  CORUM                  Gata1-Snf2h complex
    ## 6      0.02 1.0000000 CORUM:5547  CORUM                   Cdc2-Ccnb1 complex
    ##   effective_domain_size source_order       parents
    ## 1                  1073          332 CORUM:0000000
    ## 2                  1073           90 CORUM:0000000
    ## 3                  1073          330 CORUM:0000000
    ## 4                  1073          250 CORUM:0000000
    ## 5                  1073          236 CORUM:0000000
    ## 6                  1073          437 CORUM:0000000

The result data.frame contains the following columns:

    names(gostres$meta)

    ## [1] "query_metadata"  "result_metadata" "genes_metadata"  "timestamp"      
    ## [5] "version"

## Multiple queries

The function `gost()` also allows to perform enrichment on multiple
input gene lists. Multiple queries are automatically detected if the
input query is a list of vectors with gene identifiers and the results
are combined into identical data.frame as in case of single query.

    multi_gostres1 <- gost(query = list("chromX" = c("X:1000:1000000", "rs17396340", 
                                                     "GO:0005005", "ENSG00000156103", "NLRP1"),
                                 "chromY" = c("Y:1:10000000", "rs17396340", 
                                              "GO:0005005", "ENSG00000156103", "NLRP1")), 
                           multi_query = FALSE)

    head(multi_gostres1$result, 3)

    ##    query significant      p_value term_size query_size intersection_size
    ## 1 chromX        TRUE 2.737754e-30        88         22                16
    ## 2 chromX        TRUE 9.341439e-22       283         22                16
    ## 3 chromX        TRUE 9.896007e-22       284         22                16
    ##   precision     recall    term_id source                         term_name
    ## 1 0.7272727 0.18181818 GO:0048013  GO:BP ephrin receptor signaling pathway
    ## 2 0.7272727 0.05653710 GO:0007411  GO:BP                     axon guidance
    ## 3 0.7272727 0.05633803 GO:0097485  GO:BP        neuron projection guidance
    ##   effective_domain_size source_order                            parents
    ## 1                 18123        14536                         GO:0007169
    ## 2                 18123         3298             GO:0007409, GO:0097485
    ## 3                 18123        21913 GO:0006928, GO:0006935, GO:0048812

The column “query” in the result dataframe will now contain the
corresponding name for the query. If no name is specified, then the
query name is defined as the order of query with the prefix “query\_.”
Another option for multiple gene lists is setting the parameter
`multiquery = TRUE`. Then the results from all of the input queries are
grouped according to term IDs for better comparison.

    multi_gostres2 <- gost(query = list("chromX" = c("X:1000:1000000", "rs17396340",
                                                     "GO:0005005", "ENSG00000156103", "NLRP1"),
                                 "chromY" = c("Y:1:10000000", "rs17396340", 
                                              "GO:0005005", "ENSG00000156103", "NLRP1")), 
                           multi_query = TRUE)

    head(multi_gostres2$result, 3)

    ##              term_id                   p_values significant term_size
    ## 1         GO:0005005 7.123326e-48, 2.269088e-44  TRUE, TRUE        16
    ## 2         GO:0005003 6.895195e-45, 2.193428e-41  TRUE, TRUE        19
    ## 3 REAC:R-HSA-3928665 3.712329e-33, 1.629532e-32  TRUE, TRUE        49
    ##   query_sizes intersection_sizes source                              term_name
    ## 1      23, 32             16, 16  GO:MF transmembrane-ephrin receptor activity
    ## 2      23, 32             16, 16  GO:MF               ephrin receptor activity
    ## 3      20, 21             16, 16   REAC EPH-ephrin mediated repulsion of cells
    ##   effective_domain_size source_order            parents
    ## 1                 18679         1536         GO:0005003
    ## 2                 18679         1534         GO:0004714
    ## 3                 10622          753 REAC:R-HSA-2682334

## Visualization

The enrichment results are visualized with a Manhattan-like-plot using
the function `gostplot()` and the previously found gost results
*gostres*:

    #gostplot(gostres, capped = TRUE, interactive = TRUE)

The function `publish_gostplot()` takes the static plot object as an
input and enables to highlight a selection of interesting terms from the
results with numbers and table of results. These can be set with
parameter highlight\_terms listing the term IDs in a vector or as a
data.frame with column “term\_id” such as a subset of the result
dataframe.

First we create the static plot

    p <- gostplot(gostres, capped = FALSE, interactive = FALSE)
    p

<img src="scRNAseq_Seurat_files/figure-markdown_strict/gost_static-1.png" style="display: block; margin: auto;" />

Then we make it in high quality. We can add highlighted terms if we want
with the *highlight\_terms* argument.

    pp <- publish_gostplot(p, highlight_terms = c("CORUM:3047"),
                           width = NA, height = NA, filename = NULL )

<img src="scRNAseq_Seurat_files/figure-markdown_strict/gost_publishing-1.png" style="display: block; margin: auto;" />

The gost results can also be visualized with a table. The
`publish_gosttable()` function will create a nice-looking table with the
result statistics for the highlight\_terms from the result data.frame.
The highlight\_terms can be a vector of term IDs or a subset of the
results.

    publish_gosttable(gostres, highlight_terms = gostres$result[c(1:10),],
                            use_colors = TRUE, 
                            show_columns = c("source", "term_name", "term_size", "intersection_size"),
                            filename = NULL) 

<img src="scRNAseq_Seurat_files/figure-markdown_strict/gost_table-1.png" style="display: block; margin: auto;" />

# Pseudotime

We will use the package
[Slignshot](https://bioconductor.org/packages/release/bioc/vignettes/slingshot/inst/doc/vignette.html).
Slingshot was designed to model developmental trajectories in
single-cell RNA sequencing data and serve as a component in an analysis
pipeline after dimensionality reduction and clustering.

    adata_SCE <- as.SingleCellExperiment(adata)
    sde <- slingshot(adata_SCE, reducedDim = 'PCA')

The using the `slingshot` package it is possible to plot our components
and the calculated trajectory. Unfortunately, it does not work with
ggplot.

    cols <- colorRampPalette(viridis::viridis(3))
    cols <- cols(50)[cut(sde$slingPseudotime_1, breaks=50)]
    plot(reducedDim(adata_SCE, "PCA"), pch=16, col = cols) # adata_SCE contains, PC1 and PC2, and colored by pseudotime
    lines(slingshot::SlingshotDataSet(sde), lwd = 2) # here we plot the trajectory

<img src="scRNAseq_Seurat_files/figure-markdown_strict/sling_pseudo_plot-1.png" style="display: block; margin: auto;" />

It is possible to create a `FeaturePlot()` with both pseudotime and the
trajectory by extracting the information from our `slingshot` object and
add it to our `Seurat` object, but it is a bit inconvenient. First,
pseudotime can be extracted with the function `slingPseudotime()`. The
trajectory curves are stored inside the slingshot object `sde`; in this
case, it is easier since there is only one trajectory.

    adata$pseudotime <- slingPseudotime(sde)
    trajectory <- data.frame(SlingshotDataSet(sde)@curves[["Lineage1"]][["s"]])
    head(trajectory)

    ##       PC_1       PC_2     PC_3      PC_4      PC_5       PC_6       PC_7
    ## 1 16.35562 -11.962698 16.15586 -4.601536 -5.039318 0.31009288 -1.0479768
    ## 2 16.04464 -11.275443 15.02745 -4.182979 -4.473788 0.25979631 -0.9271683
    ## 3 15.73366 -10.588173 13.89902 -3.764424 -3.908291 0.20950282 -0.8063613
    ## 4 15.42250  -9.900414 12.77014 -3.345910 -3.343873 0.15930534 -0.6856004
    ## 5 15.11072  -9.211110 11.63987 -2.927527 -2.782916 0.10932772 -0.5650254
    ## 6 14.79767  -8.518681 10.50682 -2.509419 -2.229046 0.05966221 -0.4449966
    ##        PC_8        PC_9      PC_10     PC_11      PC_12      PC_13       PC_14
    ## 1 0.3254987 -0.07134998 -0.7871131 0.8965533 -0.7646160 0.11605931 -0.05243540
    ## 2 0.2884276 -0.04177932 -0.6936638 0.7770207 -0.6634016 0.10375284 -0.04891362
    ## 3 0.2513555 -0.01221806 -0.6002311 0.6575134 -0.5622026 0.09144860 -0.04539230
    ## 4 0.2142538  0.01703939 -0.5073295 0.5388099 -0.4615034 0.07921579 -0.04188852
    ## 5 0.1770835  0.04533347 -0.4160275 0.4224308 -0.3623267 0.06719155 -0.03849633
    ## 6 0.1398380  0.07156764 -0.3281833 0.3106164 -0.2661960 0.05559170 -0.03536768
    ##       PC_15       PC_16     PC_17     PC_18      PC_19       PC_20      PC_21
    ## 1 0.2934881 -0.17907594 0.4235500 0.3886445 0.24248844 -0.10529320 0.08863358
    ## 2 0.2616397 -0.16045868 0.3663639 0.3461969 0.20809433 -0.09242672 0.07232385
    ## 3 0.2297896 -0.14183801 0.3091900 0.3037561 0.17371050 -0.07957034 0.05602695
    ## 4 0.1978901 -0.12311010 0.2524158 0.2615234 0.13965827 -0.06702953 0.04013605
    ## 5 0.1660107 -0.10410463 0.1969324 0.2197808 0.10662611 -0.05527745 0.02534756
    ## 6 0.1341776 -0.08466072 0.1440259 0.1789958 0.07558997 -0.04492861 0.01263065
    ##        PC_22      PC_23      PC_24      PC_25      PC_26       PC_27      PC_28
    ## 1 0.08889671 0.10737870 0.16991903 0.06071445 0.07491283 -0.01762855 0.13740571
    ## 2 0.08011988 0.09651028 0.14732173 0.05607936 0.06329914 -0.01652658 0.12212549
    ## 3 0.07134151 0.08563549 0.12472935 0.05143975 0.05168646 -0.01542325 0.10684189
    ## 4 0.06251334 0.07456311 0.10228788 0.04665799 0.04011298 -0.01427594 0.09146331
    ## 5 0.05354290 0.06303155 0.08016147 0.04151141 0.02879322 -0.01298514 0.07606155
    ## 6 0.04443393 0.05070978 0.05854566 0.03574929 0.01818949 -0.01142445 0.06076868
    ##        PC_29       PC_30       PC_31       PC_32       PC_33      PC_34
    ## 1 0.12421169 -0.12040235 -0.08514957 -0.01833432 -0.05863020 0.15169391
    ## 2 0.11089818 -0.10387949 -0.07543986 -0.02065316 -0.05396725 0.13211851
    ## 3 0.09758189 -0.08735898 -0.06572574 -0.02296483 -0.04930408 0.11254913
    ## 4 0.08418526 -0.07091984 -0.05588027 -0.02503999 -0.04462607 0.09316837
    ## 5 0.07072997 -0.05483307 -0.04583502 -0.02628717 -0.03975497 0.07427303
    ## 6 0.05738486 -0.03946295 -0.03560769 -0.02593073 -0.03441431 0.05630623
    ##        PC_35      PC_36       PC_37      PC_38      PC_39       PC_40
    ## 1 0.10531852 0.03937051 -0.16130429 0.01370460 0.06984152 -0.09911683
    ## 2 0.09218033 0.03805641 -0.14279212 0.01281096 0.06308274 -0.08761607
    ## 3 0.07903989 0.03673776 -0.12427754 0.01192263 0.05632439 -0.07611601
    ## 4 0.06584141 0.03526996 -0.10568779 0.01119295 0.04958260 -0.06462963
    ## 5 0.05273105 0.03329917 -0.08693079 0.01070667 0.04293414 -0.05301380
    ## 6 0.03981979 0.03034186 -0.06815030 0.01040570 0.03641309 -0.04107223
    ##         PC_41         PC_42      PC_43       PC_44       PC_45      PC_46
    ## 1 -0.06242096  0.0085915291 0.12679393 -0.09390764 -0.21561205 0.01939662
    ## 2 -0.05329910  0.0069205683 0.11126941 -0.08400701 -0.18858069 0.01822608
    ## 3 -0.04419030  0.0052483723 0.09574446 -0.07411190 -0.16155028 0.01705700
    ## 4 -0.03549195  0.0035376232 0.08019930 -0.06438479 -0.13455048 0.01593352
    ## 5 -0.02783877  0.0017374348 0.06446992 -0.05497856 -0.10765005 0.01492262
    ## 6 -0.02188553 -0.0001414273 0.04841770 -0.04606122 -0.08111784 0.01407998
    ##         PC_47         PC_48       PC_49       PC_50
    ## 1 -0.15191992 -0.0045559770 -0.04749615 -0.10971005
    ## 2 -0.13000420 -0.0031502924 -0.04274298 -0.09485364
    ## 3 -0.10810006 -0.0017433838 -0.03799244 -0.07999711
    ## 4 -0.08656298 -0.0003050211 -0.03332203 -0.06513857
    ## 5 -0.06603371  0.0010754775 -0.02881896 -0.05029519
    ## 6 -0.04733289  0.0022401686 -0.02466975 -0.03533874

Now we can use Seurat functions to plot the pseudotime and the
calculated trajectory!

    FeaturePlot(adata, features = "pseudotime", reduction = "pca") + scale_color_viridis_c() + geom_path(data = trajectory, aes(x = PC_1, y = PC_2))

<img src="scRNAseq_Seurat_files/figure-markdown_strict/pseudo_plot-1.png" style="display: block; margin: auto;" />

    DimPlot(adata, group.by = "Cell_type", reduction = "pca")

<img src="scRNAseq_Seurat_files/figure-markdown_strict/unnamed-chunk-32-1.png" style="display: block; margin: auto;" />

We can create a plot that sorts our cells by pseudotime and check that
it corresponds to their cell type:

    adata_SCE$pseudotime <- as.numeric(slingPseudotime(sde))
    ggplot(as.data.frame(colData(adata_SCE)), aes(x = pseudotime,
                                                 y = Cell_type,
                                                 colour = Cell_type)) +
      ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
      theme_classic() +
      xlab("Slingshot pseudotime") + ylab("Timepoint") +
      ggtitle("Cells ordered by Slingshot pseudotime")

<img src="scRNAseq_Seurat_files/figure-markdown_strict/pseudo_order-1.png" style="display: block; margin: auto;" />

# Session info

Finally, we create a `session_info()` table that will allow anyone to
check what versions of R and packages are we using for reproducibility
purposes.

    devtools::session_info()

    ## ─ Session info ───────────────────────────────────────────────────────────────
    ##  setting  value
    ##  version  R version 4.1.0 (2021-05-18)
    ##  os       macOS Big Sur 10.16
    ##  system   x86_64, darwin17.0
    ##  ui       X11
    ##  language (EN)
    ##  collate  en_US.UTF-8
    ##  ctype    en_US.UTF-8
    ##  tz       Europe/Copenhagen
    ##  date     2022-02-02
    ##  pandoc   2.14.0.1 @ /usr/local/bin/ (via rmarkdown)
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  package              * version  date (UTC) lib source
    ##  abind                  1.4-5    2016-07-21 [1] CRAN (R 4.1.0)
    ##  AnnotationDbi          1.56.2   2021-11-09 [1] Bioconductor
    ##  assertthat             0.2.1    2019-03-21 [1] CRAN (R 4.1.0)
    ##  backports              1.4.1    2021-12-13 [1] CRAN (R 4.1.0)
    ##  beeswarm               0.4.0    2021-06-01 [1] CRAN (R 4.1.0)
    ##  Biobase              * 2.54.0   2021-10-26 [1] Bioconductor
    ##  BiocFileCache          2.2.0    2021-10-26 [1] Bioconductor
    ##  BiocGenerics         * 0.40.0   2021-10-26 [1] Bioconductor
    ##  BiocManager            1.30.16  2021-06-15 [1] CRAN (R 4.1.0)
    ##  biomaRt              * 2.50.1   2021-11-21 [1] Bioconductor
    ##  Biostrings             2.62.0   2021-10-26 [1] Bioconductor
    ##  bit                    4.0.4    2020-08-04 [1] CRAN (R 4.1.0)
    ##  bit64                  4.0.5    2020-08-30 [1] CRAN (R 4.1.0)
    ##  bitops                 1.0-7    2021-04-24 [1] CRAN (R 4.1.0)
    ##  blob                   1.2.2    2021-07-23 [1] CRAN (R 4.1.0)
    ##  broom                  0.7.10   2021-10-31 [1] CRAN (R 4.1.0)
    ##  cachem                 1.0.6    2021-08-19 [1] CRAN (R 4.1.0)
    ##  callr                  3.7.0    2021-04-20 [1] CRAN (R 4.1.0)
    ##  cellranger             1.1.0    2016-07-27 [1] CRAN (R 4.1.0)
    ##  cli                    3.1.0    2021-10-27 [1] CRAN (R 4.1.0)
    ##  cluster                2.1.2    2021-04-17 [1] CRAN (R 4.1.0)
    ##  codetools              0.2-18   2020-11-04 [1] CRAN (R 4.1.0)
    ##  colorspace             2.0-2    2021-06-24 [1] CRAN (R 4.1.0)
    ##  cowplot                1.1.1    2020-12-30 [1] CRAN (R 4.1.0)
    ##  crayon                 1.4.2    2021-10-29 [1] CRAN (R 4.1.0)
    ##  curl                   4.3.2    2021-06-23 [1] CRAN (R 4.1.0)
    ##  data.table             1.14.2   2021-09-27 [1] CRAN (R 4.1.0)
    ##  DBI                    1.1.1    2021-01-15 [1] CRAN (R 4.1.0)
    ##  dbplyr                 2.1.1    2021-04-06 [1] CRAN (R 4.1.0)
    ##  DelayedArray           0.20.0   2021-10-26 [1] Bioconductor
    ##  deldir                 1.0-6    2021-10-23 [1] CRAN (R 4.1.0)
    ##  desc                   1.4.0    2021-09-28 [1] CRAN (R 4.1.0)
    ##  devtools               2.4.3    2021-11-30 [1] CRAN (R 4.1.0)
    ##  digest                 0.6.29   2021-12-01 [1] CRAN (R 4.1.0)
    ##  dplyr                * 1.0.7    2021-06-18 [1] CRAN (R 4.1.0)
    ##  ellipsis               0.3.2    2021-04-29 [1] CRAN (R 4.1.0)
    ##  evaluate               0.14     2019-05-28 [1] CRAN (R 4.1.0)
    ##  fansi                  0.5.0    2021-05-25 [1] CRAN (R 4.1.0)
    ##  farver                 2.1.0    2021-02-28 [1] CRAN (R 4.1.0)
    ##  fastmap                1.1.0    2021-01-25 [1] CRAN (R 4.1.0)
    ##  filelock               1.0.2    2018-10-05 [1] CRAN (R 4.1.0)
    ##  fitdistrplus           1.1-6    2021-09-28 [1] CRAN (R 4.1.0)
    ##  forcats              * 0.5.1    2021-01-27 [1] CRAN (R 4.1.0)
    ##  fs                     1.5.2    2021-12-08 [1] CRAN (R 4.1.0)
    ##  future                 1.23.0   2021-10-31 [1] CRAN (R 4.1.0)
    ##  future.apply           1.8.1    2021-08-10 [1] CRAN (R 4.1.0)
    ##  generics               0.1.1    2021-10-25 [1] CRAN (R 4.1.0)
    ##  GenomeInfoDb         * 1.30.0   2021-10-26 [1] Bioconductor
    ##  GenomeInfoDbData       1.2.7    2021-11-16 [1] Bioconductor
    ##  GenomicRanges        * 1.46.1   2021-11-18 [1] Bioconductor
    ##  ggbeeswarm             0.6.0    2017-08-07 [1] CRAN (R 4.1.0)
    ##  ggplot2              * 3.3.5    2021-06-25 [1] CRAN (R 4.1.0)
    ##  ggrepel                0.9.1    2021-01-15 [1] CRAN (R 4.1.0)
    ##  ggridges               0.5.3    2021-01-08 [1] CRAN (R 4.1.0)
    ##  globals                0.14.0   2020-11-22 [1] CRAN (R 4.1.0)
    ##  glue                   1.5.1    2021-11-30 [1] CRAN (R 4.1.0)
    ##  goftest                1.2-3    2021-10-07 [1] CRAN (R 4.1.0)
    ##  gprofiler2           * 0.2.1    2021-08-23 [1] CRAN (R 4.1.0)
    ##  gridExtra              2.3      2017-09-09 [1] CRAN (R 4.1.0)
    ##  gtable                 0.3.0    2019-03-25 [1] CRAN (R 4.1.0)
    ##  haven                  2.4.3    2021-08-04 [1] CRAN (R 4.1.0)
    ##  highr                  0.9      2021-04-16 [1] CRAN (R 4.1.0)
    ##  hms                    1.1.1    2021-09-26 [1] CRAN (R 4.1.0)
    ##  htmltools              0.5.2    2021-08-25 [1] CRAN (R 4.1.0)
    ##  htmlwidgets            1.5.4    2021-09-08 [1] CRAN (R 4.1.0)
    ##  httpuv                 1.6.3    2021-09-09 [1] CRAN (R 4.1.0)
    ##  httr                   1.4.2    2020-07-20 [1] CRAN (R 4.1.0)
    ##  ica                    1.0-2    2018-05-24 [1] CRAN (R 4.1.0)
    ##  igraph                 1.2.9    2021-11-23 [1] CRAN (R 4.1.0)
    ##  IRanges              * 2.28.0   2021-10-26 [1] Bioconductor
    ##  irlba                  2.3.5    2021-12-06 [1] CRAN (R 4.1.0)
    ##  jsonlite               1.7.2    2020-12-09 [1] CRAN (R 4.1.0)
    ##  KEGGREST               1.34.0   2021-10-26 [1] Bioconductor
    ##  KernSmooth             2.23-20  2021-05-03 [1] CRAN (R 4.1.0)
    ##  knitr                  1.36     2021-09-29 [1] CRAN (R 4.1.0)
    ##  labeling               0.4.2    2020-10-20 [1] CRAN (R 4.1.0)
    ##  later                  1.3.0    2021-08-18 [1] CRAN (R 4.1.0)
    ##  lattice                0.20-45  2021-09-22 [1] CRAN (R 4.1.0)
    ##  lazyeval               0.2.2    2019-03-15 [1] CRAN (R 4.1.0)
    ##  leiden                 0.3.9    2021-07-27 [1] CRAN (R 4.1.0)
    ##  lifecycle              1.0.1    2021-09-24 [1] CRAN (R 4.1.0)
    ##  limma                  3.50.0   2021-10-26 [1] Bioconductor
    ##  listenv                0.8.0    2019-12-05 [1] CRAN (R 4.1.0)
    ##  lmtest                 0.9-39   2021-11-07 [1] CRAN (R 4.1.0)
    ##  lubridate              1.8.0    2021-10-07 [1] CRAN (R 4.1.0)
    ##  magrittr               2.0.1    2020-11-17 [1] CRAN (R 4.1.0)
    ##  MASS                   7.3-54   2021-05-03 [1] CRAN (R 4.1.0)
    ##  Matrix                 1.4-0    2021-12-08 [1] CRAN (R 4.1.0)
    ##  MatrixGenerics       * 1.6.0    2021-10-26 [1] Bioconductor
    ##  matrixStats          * 0.61.0   2021-09-17 [1] CRAN (R 4.1.0)
    ##  memoise                2.0.1    2021-11-26 [1] CRAN (R 4.1.0)
    ##  mgcv                   1.8-38   2021-10-06 [1] CRAN (R 4.1.0)
    ##  mime                   0.12     2021-09-28 [1] CRAN (R 4.1.0)
    ##  miniUI                 0.1.1.1  2018-05-18 [1] CRAN (R 4.1.0)
    ##  modelr                 0.1.8    2020-05-19 [1] CRAN (R 4.1.0)
    ##  munsell                0.5.0    2018-06-12 [1] CRAN (R 4.1.0)
    ##  nlme                   3.1-153  2021-09-07 [1] CRAN (R 4.1.0)
    ##  parallelly             1.29.0   2021-11-21 [1] CRAN (R 4.1.0)
    ##  patchwork            * 1.1.1    2020-12-17 [1] CRAN (R 4.1.0)
    ##  pbapply                1.5-0    2021-09-16 [1] CRAN (R 4.1.0)
    ##  pillar                 1.6.4    2021-10-18 [1] CRAN (R 4.1.0)
    ##  pkgbuild               1.3.0    2021-12-09 [1] CRAN (R 4.1.0)
    ##  pkgconfig              2.0.3    2019-09-22 [1] CRAN (R 4.1.0)
    ##  pkgload                1.2.4    2021-11-30 [1] CRAN (R 4.1.0)
    ##  plotly                 4.10.0   2021-10-09 [1] CRAN (R 4.1.0)
    ##  plyr                   1.8.6    2020-03-03 [1] CRAN (R 4.1.0)
    ##  png                    0.1-7    2013-12-03 [1] CRAN (R 4.1.0)
    ##  polyclip               1.10-0   2019-03-14 [1] CRAN (R 4.1.0)
    ##  prettyunits            1.1.1    2020-01-24 [1] CRAN (R 4.1.0)
    ##  princurve            * 2.1.6    2021-01-18 [1] CRAN (R 4.1.0)
    ##  processx               3.5.2    2021-04-30 [1] CRAN (R 4.1.0)
    ##  progress               1.2.2    2019-05-16 [1] CRAN (R 4.1.0)
    ##  promises               1.2.0.1  2021-02-11 [1] CRAN (R 4.1.0)
    ##  ps                     1.6.0    2021-02-28 [1] CRAN (R 4.1.0)
    ##  purrr                * 0.3.4    2020-04-17 [1] CRAN (R 4.1.0)
    ##  R.methodsS3            1.8.1    2020-08-26 [1] CRAN (R 4.1.0)
    ##  R.oo                   1.24.0   2020-08-26 [1] CRAN (R 4.1.0)
    ##  R.utils                2.11.0   2021-09-26 [1] CRAN (R 4.1.0)
    ##  R6                     2.5.1    2021-08-19 [1] CRAN (R 4.1.0)
    ##  RANN                   2.6.1    2019-01-08 [1] CRAN (R 4.1.0)
    ##  rappdirs               0.3.3    2021-01-31 [1] CRAN (R 4.1.0)
    ##  RColorBrewer           1.1-2    2014-12-07 [1] CRAN (R 4.1.0)
    ##  Rcpp                   1.0.7    2021-07-07 [1] CRAN (R 4.1.0)
    ##  RcppAnnoy              0.0.19   2021-07-30 [1] CRAN (R 4.1.0)
    ##  RCurl                  1.98-1.5 2021-09-17 [1] CRAN (R 4.1.0)
    ##  readr                * 2.1.1    2021-11-30 [1] CRAN (R 4.1.0)
    ##  readxl                 1.3.1    2019-03-13 [1] CRAN (R 4.1.0)
    ##  remotes                2.4.2    2021-11-30 [1] CRAN (R 4.1.0)
    ##  reprex                 2.0.1    2021-08-05 [1] CRAN (R 4.1.0)
    ##  reshape2               1.4.4    2020-04-09 [1] CRAN (R 4.1.0)
    ##  reticulate             1.22     2021-09-17 [1] CRAN (R 4.1.0)
    ##  rhdf5                * 2.38.0   2021-10-26 [1] Bioconductor
    ##  rhdf5filters           1.6.0    2021-10-26 [1] Bioconductor
    ##  Rhdf5lib               1.16.0   2021-10-26 [1] Bioconductor
    ##  rlang                  0.4.12   2021-10-18 [1] CRAN (R 4.1.0)
    ##  rmarkdown              2.11     2021-09-14 [1] CRAN (R 4.1.0)
    ##  ROCR                   1.0-11   2020-05-02 [1] CRAN (R 4.1.0)
    ##  rpart                  4.1-15   2019-04-12 [1] CRAN (R 4.1.0)
    ##  rprojroot              2.0.2    2020-11-15 [1] CRAN (R 4.1.0)
    ##  RSpectra               0.16-0   2019-12-01 [1] CRAN (R 4.1.0)
    ##  RSQLite                2.2.9    2021-12-06 [1] CRAN (R 4.1.0)
    ##  rstudioapi             0.13     2020-11-12 [1] CRAN (R 4.1.0)
    ##  rsvd                   1.0.5    2021-04-16 [1] CRAN (R 4.1.0)
    ##  Rtsne                  0.15     2018-11-10 [1] CRAN (R 4.1.0)
    ##  rvest                  1.0.2    2021-10-16 [1] CRAN (R 4.1.0)
    ##  S4Vectors            * 0.32.3   2021-11-21 [1] Bioconductor
    ##  scales                 1.1.1    2020-05-11 [1] CRAN (R 4.1.0)
    ##  scattermore            0.7      2020-11-24 [1] CRAN (R 4.1.0)
    ##  sctransform            0.3.2    2020-12-16 [1] CRAN (R 4.1.0)
    ##  sessioninfo            1.2.2    2021-12-06 [1] CRAN (R 4.1.0)
    ##  Seurat               * 4.0.5    2021-10-17 [1] CRAN (R 4.1.0)
    ##  SeuratObject         * 4.0.4    2021-11-23 [1] CRAN (R 4.1.0)
    ##  SeuratWrappers       * 0.3.0    2021-06-22 [1] Github (satijalab/seurat-wrappers@fba0662)
    ##  shiny                  1.7.1    2021-10-02 [1] CRAN (R 4.1.0)
    ##  SingleCellExperiment * 1.16.0   2021-10-26 [1] Bioconductor
    ##  slingshot            * 2.2.0    2021-10-26 [1] Bioconductor
    ##  spatstat.core          2.3-2    2021-11-26 [1] CRAN (R 4.1.0)
    ##  spatstat.data          2.1-0    2021-03-21 [1] CRAN (R 4.1.0)
    ##  spatstat.geom          2.3-1    2021-12-10 [1] CRAN (R 4.1.0)
    ##  spatstat.sparse        2.0-0    2021-03-16 [1] CRAN (R 4.1.0)
    ##  spatstat.utils         2.3-0    2021-12-12 [1] CRAN (R 4.1.0)
    ##  stringi                1.7.6    2021-11-29 [1] CRAN (R 4.1.0)
    ##  stringr              * 1.4.0    2019-02-10 [1] CRAN (R 4.1.0)
    ##  SummarizedExperiment * 1.24.0   2021-10-26 [1] Bioconductor
    ##  survival               3.2-13   2021-08-24 [1] CRAN (R 4.1.0)
    ##  tensor                 1.5      2012-05-05 [1] CRAN (R 4.1.0)
    ##  testthat               3.1.1    2021-12-03 [1] CRAN (R 4.1.0)
    ##  tibble               * 3.1.6    2021-11-07 [1] CRAN (R 4.1.0)
    ##  tidyr                * 1.1.4    2021-09-27 [1] CRAN (R 4.1.0)
    ##  tidyselect             1.1.1    2021-04-30 [1] CRAN (R 4.1.0)
    ##  tidyverse            * 1.3.1    2021-04-15 [1] CRAN (R 4.1.0)
    ##  TrajectoryUtils      * 1.2.0    2021-10-26 [1] Bioconductor
    ##  tzdb                   0.2.0    2021-10-27 [1] CRAN (R 4.1.0)
    ##  usethis                2.1.5    2021-12-09 [1] CRAN (R 4.1.0)
    ##  utf8                   1.2.2    2021-07-24 [1] CRAN (R 4.1.0)
    ##  uwot                   0.1.11   2021-12-02 [1] CRAN (R 4.1.0)
    ##  vctrs                  0.3.8    2021-04-29 [1] CRAN (R 4.1.0)
    ##  vipor                  0.4.5    2017-03-22 [1] CRAN (R 4.1.0)
    ##  viridis                0.6.2    2021-10-13 [1] CRAN (R 4.1.0)
    ##  viridisLite            0.4.0    2021-04-13 [1] CRAN (R 4.1.0)
    ##  withr                  2.4.3    2021-11-30 [1] CRAN (R 4.1.0)
    ##  xfun                   0.28     2021-11-04 [1] CRAN (R 4.1.0)
    ##  XML                    3.99-0.8 2021-09-17 [1] CRAN (R 4.1.0)
    ##  xml2                   1.3.3    2021-11-30 [1] CRAN (R 4.1.0)
    ##  xtable                 1.8-4    2019-04-21 [1] CRAN (R 4.1.0)
    ##  XVector                0.34.0   2021-10-26 [1] Bioconductor
    ##  yaml                   2.2.1    2020-02-01 [1] CRAN (R 4.1.0)
    ##  zlibbioc               1.40.0   2021-10-26 [1] Bioconductor
    ##  zoo                    1.8-9    2021-03-09 [1] CRAN (R 4.1.0)
    ## 
    ##  [1] /Library/Frameworks/R.framework/Versions/4.1/Resources/library
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
