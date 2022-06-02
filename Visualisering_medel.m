%% Visualusering medelvärde
clear all; close all; clc;


%%
rate_ana_1 = [0
0
0
0.006201550388
0.008965517241
0.009971004201
0.01060979132
0.009156036736
0.01065946904
0.01250275264
0.01232991877
0.009540423889
0.01420952671
0.01466703534
0.02117264507
0.02210440971
0.01242017484
0.02384740547
0.02697743603
0.03385213769
0.02723128287
0.02168340352
0.02127397324
0.02834079773
0.02963485859
0.03227028954
0.03092378504
0.02389149825
0.02312930253
0.03244218871
0.02337547247
0.02544462435
0.02523979114
0.0124894239
0.01732162776
0.0131725433
0.01535553139
0.01648146779
0.01472020856
0.01668394562
0.01167866751
0.01030656071
0.02598889919
0.03730893782
0.03325969993
0.01585850262
0.0281700808
0.02001445271
0.01961380269
0.02262699943
0.02073424671
0.02189573404
0.02378131587
0.02160652764
0.01680788698
0.01552767661
0.02671365392
0.01999926268
0.02619192203
0.01702947076
0.01957844813
0.02083741087
0.02122110722
0.015331817
0.01342245786
0.01936170286
0.0230240955
0.02881865093
0.03627110823
0.02640251529
0.02098179237
0.0150285872
0.01166795642
0.0110706695
0.009919549083
0.01051426889
0.01389136131
0.0215035066
0.01959536167
0.02930688065
0.02643173012
0.02123636012
0.01687179128
0.02939380714
0.03236402541
0.03631713469
0.03095269874
0.03995778411
0.0296465683
0.03328313704
0.03552008421
0.03400237402
0.03040696774
0.01371716258
0.01968929599
0.01896639104
0.02180486471
0.03459675477
0.03175781873
0.03098736404
0.03378564333
0.02847600335
0.03203129307
0.03761910325
0.03229796542
0.04114467298
0.0318715929
0.03245240491
0.03478619724
0.0424181488
0.03808687306
0.04003271783
0.03781857409
0.04007407329
0.04198738268
0.03681311784
0.04073863595
0.03899914176
0.03456797971
0.03192257014
0.04525048827
0.03259525437
0.02916750644
0.02566161584
0.03463562133
0.0401590882
0.02915503821
0.03047015352
0.03591870788
0.04148001467
0.04632723575
0.04643390115
0.0405745799
0.05502556191
0.05950049505
0.05108337327
0.05226543541
0.05768815856
0.05768954426
0.05866232021
0.05816553803
0.05476930153
0.05731458305
0.05932818208
0.0601741798
0.05306028582
0.0479834903
0.04249582561
0.03819250142
0.03936162121
0.04528010715
0.04970351532
0.04412762122
0.04132379014
0.05019078965
0.04924061191
0.04849657185
0.04140140878
0.04681009039
0.04849558209
0.05660966953
0.06020357293
0.06651363905
0.06436698636
0.06370063672
0.05530185313
0.06163260989
0.06186255919
0.05855143316
0.05607696899
0.06185228205
0.06260689983
0.06683561395
0.07232314586
0.08482116987
0.08186596746
0.09386691574
0.07892687884
0.0983795455
0.09946903999
0.112173936
0.1114253882
0.1197180851
0.1145460621
0.110296187
0.09823725927
0.09378469502
0.0774697313
0.08167388602
0.07876604234
0.07978242157
0.07936360993
0.07698661386
0.07701763411
0.07601491167
0.07699020327
0.0790165773
0.081381399
0.08057280639
0.08599305872
0.08665069889
0.08745996681
0.09287339665
0.09112816965
0.08642161332
0.08317915598
0.08451481087
0.08732746913
0.08500234846
0.08210231789
0.07369070013
0.07778896089
0.08495352722
0.08357862221
0.08564350565
0.09001392256
0.09845020561
0.1006946732
0.1026576589
0.1015095214
0.08838466706
0.08867282191
0.09005903603
0.08556528231
0.09365216551
0.09424130338
0.0986668796
0.09640665273
0.1051798442
0.0994579828
0.09695307573
0.09971830052
0.09952130198
0.09839310944
0.1020985381
0.09802621941
0.09246513245
0.09325435372
0.09739489814
0.09780388217
0.0887515171
0.08568687608
0.0908321212
0.08689623493
0.0855824388
0.08707746537
0.08899981032
0.08565212475
0.08655985163
0.08288736894
0.08384069411
0.08470497525
0.08922198361
0.09074364094
0.08846325404
0.08681688658
0.08616994202
0.08566951102
0.08966900283
0.08532069279
0.08738138925
0.08708929223
0.08705783803
0.08593884511
0.08288072787
0.08145322863
0.07824482811
0.07607302564
0.07519632297
0.07377129335
0.07797512969
0.07682349434
0.08036970271
0.08362179614
0.08120947613
0.08134409046
0.0822325154
0.08355081268
0.0831215533
0.08303982832
0.08320305713
0.08036749136
0.0860107448
0.09206925904
0.08890688112
0.0893808734
0.08862750936
0.08860111451
0.08719218764
0.08212993232
0.08381140178
0.08362681951
0.08242168395
0.08099173571
0.08038227763
0.0803442166
0.07914914477
0.07736571043
0.08035156879
0.07602941042
0.07934373344
0.08050031684
0.07885871574
0.07883318736
0.07866479636
0.07843942674
0.07960993326
0.07917538221
0.07990549325
0.07962712824
0.07785893949
0.08125316065
0.08095501006
0.08001847201
0.07888716548
0.07884495892
0.07882190968
0.08036407185
0.07931872815
0.0770812292
0.07788527003
0.07458487652
0.0751295315
0.07624366871
0.07568420108
0.07767744944
0.07431338867
0.07401446595
0.07342715644
0.07266482757
0.07316358114
0.07386872143
0.07492511821
0.07483780366
0.07300134629
0.07218439343
0.07169018865
0.07212489281
0.07285148119
0.07370808076
0.07331541822
0.07297938668
0.07348291349
0.07452305956
0.07650495812
0.07661508838
0.07615888988
0.07752263802
0.07650520155
0.078449507
0.07964803053
0.0776762999
0.07691263873
0.07600424391
0.07658063633
0.07967837877
0.07944615966
0.08073968668
0.08020256209
0.08010610398
0.07598107163
0.07796172979
0.07818562922
0.07730019001
0.07704246529
0.07762233571
0.0784274008
0.07799720176
0.07807230076
0.07723691157
0.07707724641
0.07733236326
0.07781922965
0.07816897873
0.08067485454
0.07924474175
0.07949345357
0.07852962652
0.07998786001
0.08051373152
0.0817369285
0.08084277457
0.07936353509
0.07809674438
0.07619976698
0.07592863256
0.07515785432
0.07456593403
0.07380216127
0.07404896319
0.07420268945
0.07377266081
0.07528126262
0.07470044398
0.07372376991
0.07337669569
0.07453525518
0.07411725823
0.07518274791
0.07585832903
0.07660781749
0.07550957727
0.07457357676
0.07446900862
0.07365571669
0.07519766211
0.07461195126
0.07486455931
0.07371176558
0.07358805006
0.07303420029
0.07375263535
0.07388162581
0.07427343495
0.07298615014
0.0728698575
0.0717169279
0.07201426963
0.07201544672
0.07254146055
0.07254379832
0.07321062777
0.07252528868
0.07283095773
0.07172931941
0.07065516064
0.07006792971
0.07083449029
0.07184628439
0.07058301438
0.06995761307
0.07068259031
0.07174027578
0.07196774886
0.07147323999
0.07126548718
0.07128242827
0.07174966073
0.07189563459
0.07259119827
0.07290948734
0.0729255138
0.07247077253
0.07352883498
0.0733188727
0.07266936547
0.07366055625
0.03696137006
0.03681463944
0.03539028841
0.03497182188
0.03498914811
0.03460031383
0.03489371788
0.03514203881
0.04602831182
0.04581827274
0.04508780651
0.04255085678
0.03701902767
0.03718032825
0.03793575124
0.03787023121
0.03766143863
0.03762141181
0.04749080183
0.03973669726
0.03199941029
0.02551045341
0.02543946894
0.02402946321
0.02445657585
0.02444053286
0.02522747539
0.02593592577
0.03164345058
0.03419489891
0.03285969676
0.02949242954
0.02841917265
0.02784251692
0.02730293156
0.02783357836
0.02761311003
0.02774698176
0.02882423555
0.03123765977
0.02768786637
0.02116111651
0.02058496816
0.02054293163
0.0206070885
0.02129882995
0.02139955972
0.0211183186
0.02219482332
0.02484326952
0.02444543542
0.0256234388
0.02342089234
0.02299301385
0.02376047169
0.0226755139
0.02410254507
0.02290911802
0.02222905541
0.02366691581
0.02467425048
0.02393915593
0.0219168225
0.0212685886
0.02111429867
0.02099466073
0.02023010762
0.01983034631
0.02018522145
0.02075710101
0.02143924135
0.01971125686
0.02006891827
0.01969573223
0.01854968003
0.01796333594
0.01759387044
0.01776616308
0.01821319623
0.01811212723
0.01763840018
0.01778953079
0.01847083153
0.01806580194
0.01697827536
0.01668124752
0.0168141993
0.01736692364
0.01736640589
0.018006701
0.01810433466
0.01879390002
0.01828236497
0.01831778818
0.01706098196
0.01713076358
0.01720384093
0.01770643967
0.01704999455
0.01701998385
0.0178049087
0.01812434387
0.01885718276
0.01837164951
0.01785333873
0.01772763395
0.01813488801
0.01815255405
0.01816854276
0.01831444584
0.01842344729
0.01818511772
0.01810189479
0.01779610517
0.01744297789
0.01768438053
0.01777021112
0.01795019585
0.01796887523
0.01792352783
0.018018492
0.01836734321
0.01907383421
0.0185349542
0.01839622027
0.01789666131
0.01799318351
0.01834322386
0.01836354625
0.01860390058
0.01854731231
0.01856546625
0.01920491837
0.01916156507
0.01878998175
0.01835754181
0.01822236812
0.0179997327
0.01814582904
0.01874466608
0.01914239888
0.01849404068
0.01896044609
0.01921831533
0.01908293763
0.01880016225
0.01831506611
0.01787395404
0.01760657654
0.01771609587
0.01819582324
0.01808373795
0.0178967746
0.01827176227
0.01822128401
0.0183005346
0.01797908497
0.01797416741
0.01778779794
0.01820952418
0.01788426379
0.01789432294
0.01804006871
0.01852139578
0.01856831772
0.01890805001
0.01894412619
0.01936070526
0.01878454778
0.01812924631
0.01762925371
0.01766756277
0.01761363581
0.01758686262
0.01941799431
0.01969256445
0.01960069342
0.01932826778
0.01938756031
0.01934994383
0.01948170089
0.02132255622
0.02131025324
0.02125763818
0.02133347488
0.02144371531
0.02118572974
0.02076212203
0.01793420947
0.02256934722
0.02015661677
0.01806287596
0.02672556432
0.02429369292
0.02172260032
0.02582870283
0.02280787149
0.02751343586
0.02342710189
0.02747086396
0.02603703856
0.02942755327
0.0306074736
0.03559883782
0.03190506146
0.02696289676
0.03763012716
0.03980120023
0.04035613819
0.04028102875
0.04106423131
0.04206792518
0.04115500622
0.03808559014
0.03817913348
0.04111451132
0.03275557717
0.03086640374
0.03310646258
0.03305168272
0.03383526551
0.03279063734
0.03467561587
0.03062207076
0.03100095721
0.03184981902
0.03787735667
0.03709086343
0.03150271292
0.03171745888
0.03155791299
0.03165867792
0.03165972229
0.03104560992
0.03414547361
0.03211852949
0.03068140195
0.03522931479
0.03584604654
0.03090992714
0.02965581934
0.03164922342
0.03015314294
0.02711292533
0.03499211207
0.03250847078
0.03143878367
0.03418267244
0.03180639547
0.03095611058
0.03203271597
0.0330764216
0.03450822566
0.03236689062
0.0299393156
0.02915635989
0.03000272423
0.03357366987
0.02553359056
0.02979237058
0.02544972949
0.02953264175
0.02765854545
0.02347429838
0.03111405589
0.02777751946
0.03631597725
0.03469952589
0.02831901627
0.03498704847
0.03630766991
0.03131701151
0.03818469804
0.03441433263
0.043528342
0.03345488883
0.04285566714
0.04539101418
0.03620232891
0.03637255859
0.03484691611
0.0330008005
0.0284691813
0.02508776838
0.02764594523
0.03471048532
0.03969682476
0.03036534663
0.02880838988
0.02807445151
0.02709176001
0.02628765923
0.02281813016
0.02075627811
0.02071157581
0.0218059932
0.02889966405
0.02053590246
0.02309682993
0.02053903327
0.01837197428
0.02175154307
0.02024100837
0.01979993474
0.01891202667
0.01876310303
0.019285685
0.01903728426
0.01940798926
0.01870479313
0.01791000498
0.01911338741
0.01903323113
0.01784420016
0.0180270665
0.01778964254
0.01741479598
0.01979220612
0.02176834484
0.0172628964
0.01833329398
0.01790778239
0.01818214798
0.01788383949
0.0177997335
0.01994111224
0.01831962776
0.01919713354
0.01811915556
0.01731837287
0.01787776309
0.01689056387
0.01688513828
0.01711462595
0.01697876894
0.0171953451
0.01732224905
0.01763121416
0.0171368601
0.01721311761
0.01783836761
0.01774196972
0.01717922813
0.01763286196
0.01723728493
0.01775684271
0.01810537774
0.01820230652
0.02067064829
0.0178788969
0.01717245766
0.01758006838
0.01728413837
0.01755024278
0.01743225215
0.01751723558
0.01754796717
0.01784756962
0.01710764293
0.01741524055
0.01675210812
0.01665967154
0.01649984882
0.01664768233
0.01671114036
0.01786094154
0.01689675515
0.01727871046
0.0169364445
0.01701784401
0.01665130278
0.01648668964
0.01656099143
0.01656762609
0.01656219441
0.01637190361
0.01710234111
0.01748954293
0.01653643478
0.0165247952
0.01626282671
0.01634441489
0.01640601608
0.01630841915
0.01637881687
0.01684340891
0.01759526482
0.01684058785
0.01672422079
0.01646628434
0.01646246379
0.01645452647
0.01631338472
0.01662763725
0.01673612857
0.0171909931
0.01776950251
0.01756003998
0.01680488283
0.01644993909
0.01676335482
0.01677613518
0.01697210544
0.01665879147
0.01676041651
0.016775257
0.01694296389
0.01723023962
0.01685631218
0.01658868136
0.01635172371
0.01624551261
0.01629599622
0.01640274013
0.0164055369
0.01663039642
0.01694280264
0.01665271728
0.01639527462
0.01630984129
0.01634273083
0.01654019479
0.01651152648
0.01624126614
0.01626487366
0.01631104568
0.01658319855
0.0164956138
0.01650244401
0.01619345445
0.01633391285
0.01613468834
0.01634821037
0.01650421981
0.01655585856
0.01660994392
0.01675251516
0.01635155448
0.01625234127
0.01630481475
0.01630196196
0.01630303938
0.01664123877
0.01669671874
0.01663411315
0.01642431021
0.01634568705
0.01612253336
0.01615659141
0.0162632716
0.01625647286
0.01631045191
0.01640562542
0.01636984876
0.0167458009
0.01656822235
0.01694491776
0.01668278197
0.01681081034
0.01634395829
0.01646296997
0.01645766787
0.01674671773
0.01651776596
0.01663716661
0.0163313579
0.01641051054
0.01639560431
0.01673435577
0.01647193198
0.01660377035
0.01671370356
0.01683779997
0.01658768196
0.01645974115
0.01644563429
0.01642850604
0.01638817829
0.01641636642
0.01635131163
0.01620771948
0.01637633998
0.01646910013
0.01657699312
0.01673337602
0.01671484485
0.01652713701
0.01682048419
0.01654909095
0.01644705268
0.01595970355
0.01587814601
0.01621874481
0.01639373711
0.01628929077
0.01611916988
0.01620910008
0.0159155973
0.01585379415
0.01591357698
0.01593753171
0.01595358269
0.01605721667
0.01615236119
0.01611131473
0.01609159635
0.01611550638
0.01599794443
0.01600897082
0.01611255645
0.01597352268
0.01621289667
0.01629267621
0.01633618826
0.01631119571
0.01643731047
0.01629202773
0.01625266657
0.01604000215
0.0160232413
0.01594128538
0.01614962103
0.01600192701
0.01602571569
0.01601339691
0.01611142055
0.01624106222
0.01606463234
0.01591857554
0.01580508786
0.01572729242
0.01579867591
0.01614440763
0.01583886614
0.01589949451
0.01586489189
0.01591007409
0.0159129432
0.01592756436
0.01584875405
0.01607767971
0.01607902588
0.01593806533
0.01583601758
0.01599830742
0.01598833557
0.01594992239
0.0159483102
0.01602122303
0.01605383877
0.01600826258
0.01602863917
0.01585365845
0.01578342868]';
%%
rate_ana_02 = [0
0
0
0.004578754579
0.004578754579
0.01569950484
0.02128315389
0.01727494953
0.01490561486
0.01412015465
0.01232461822
0.004796456207
0.009620403385
0.008667994069
0.01197488173
0.0175179491
0.011951977
0.01220083237
0.01521121011
0.01858216851
0.02330181309
0.03024730749
0.02651804376
0.02125123211
0.01738181064
0.02409900059
0.02353380446
0.02503670979
0.01425331779
0.008394992412
0.0173403754
0.0234813537
0.02451316707
0.02029099712
0.01027858512
0.01634563077
0.02063894729
0.03428430353
0.02545962367
0.02202997351
0.02154890716
0.01849453941
0.01576136198
0.01738126263
0.02202073164
0.02575528935
0.04448141287
0.03698298513
0.04364112325
0.03899098632
0.03981352499
0.02750497303
0.023125
0.02639691125
0.02880150709
0.03883610148
0.04853616676
0.04232525125
0.04614244425
0.04000669113
0.04842415761
0.05311301634
0.0495005513
0.05242685666
0.04817960772
0.05147101047
0.058559574
0.05614521759
0.05111839294
0.04950202028
0.05442141563
0.05753599291
0.05681839499
0.06682655459
0.0590026757
0.06052773142
0.0667153896
0.06101563694
0.06302583626
0.05683839064
0.05903913512
0.0650242202
0.07324415526
0.06727748106
0.07254215555
0.07143363986
0.07344479342
0.07396125548
0.07351822691
0.07623782807
0.08701098614
0.07749089141
0.07678294312
0.07946194918
0.08364560213
0.08677256417
0.08214212821
0.08859280219
0.08938307517
0.09258153489
0.09667317351
0.1004351223
0.1033720388
0.1034923162
0.1042139252
0.1016689273
0.108423188
0.1110887202
0.1147959448
0.1233938186
0.1265169113
0.1289147704
0.1236588774
0.1153952625
0.1143915847
0.1202743068
0.1214300217
0.1154610057
0.1157870464
0.1147557564
0.1239179272
0.1269419457
0.1220487754
0.1228786865
0.1234931005
0.1147342527
0.1158249001
0.1164779086
0.1143334805
0.1171406729
0.1190953734
0.121475145
0.1189461375
0.1198113442
0.1193058661
0.1235007287
0.1218016855
0.1226754781
0.1195151865
0.1217579945
0.1239564395
0.1243643871
0.1217995803
0.1185332378
0.119242628
0.1244346235
0.1244401703
0.1262101532
0.1237709486
0.1235577395
0.1228040522
0.1277914163
0.1248177094
0.1274153614
0.1287692274
0.1318968817
0.1295055638
0.1281582469
0.1250301332
0.1237389368
0.1278452429
0.1297884354
0.1325851257
0.1287287223
0.1302201825
0.1299177682
0.1249189941
0.1262015099
0.1255053632
0.1223356773
0.1274962137
0.1276309331
0.1244157268
0.1295756613
0.1345683847
0.131462599
0.1269937584
0.1265527752
0.1266829752
0.1289032883
0.1223632509
0.1266390415
0.1280235969
0.1234428192
0.1278760458
0.1318708921
0.1291191027
0.1262381912
0.1253713728
0.1224936512
0.1232378335
0.1223391536
0.1272842503
0.1271274671
0.129664972
0.127198976
0.1305751097
0.127428673
0.1285827793
0.1332305292
0.1349740504
0.1323881961
0.1301531803
0.1317190246
0.1311644501
0.1371262332
0.1413860551
0.1342565411
0.1407783233
0.1344178405
0.1366602244
0.1349244962
0.1307020606
0.1248633916
0.1178078427
0.1236388204
0.128303113
0.130455858
0.124787445
0.1209172032
0.1246770542
0.1266898081
0.1238053904
0.1238373326
0.1233527546
0.1322127914
0.1292122257
0.1302391454
0.1230280871
0.1268862349
0.1226367676
0.1170873687
0.1182191515
0.1177002445
0.1244789246
0.1248456336
0.1184306048
0.1176652134
0.1163253802
0.1189973365
0.1138972374
0.1243454006
0.1186405423
0.1184948411
0.1239712753
0.1249263232
0.1192410642
0.1264411456
0.122692455
0.1231423888
0.1240963665
0.1262216127
0.1283234157
0.1334941774
0.1275060472
0.1305156551
0.1209172366
0.1220046824
0.1212074193
0.1198259681
0.1173983759
0.1278125433
0.1279614727
0.1378797753
0.1357589075
0.1359068276
0.1337346663
0.1341252278
0.1287147193
0.1250555182
0.1251379267
0.1317359448
0.1291007298
0.131396508
0.1311775934
0.1325927156
0.1333850002
0.1318505625
0.1300030289
0.1278132623
0.1242257355
0.1275602943
0.1252104125
0.1308523645
0.130171929
0.136645163
0.1342104511
0.1314062326
0.1329443969
0.1238782273
0.1243417421
0.1211627633
0.1287396331
0.125095634
0.1238022431
0.1257057263
0.122460676
0.1218463272
0.1214203243
0.1181383085
0.1219972843
0.1144922658
0.1140550881
0.1167786909
0.1201865419
0.1115723893
0.1121151621
0.1116711751
0.1106024097
0.1177549494
0.1072486274
0.1028199665
0.1000162741
0.1034929986
0.104316476
0.1035067963
0.1018203844
0.101727038
0.09931515884
0.09841040507
0.1004449648
0.1021339687
0.1026434262
0.1052994469
0.103010767
0.1023204114
0.101293216
0.1016750023
0.10038738
0.1037762437
0.1023268086
0.09468056781
0.09510107863
0.09906332309
0.1009594544
0.1021920806
0.09859105077
0.09403377013
0.09438638819
0.0954400169
0.09567401938
0.09177353882
0.08945228828
0.08841890106
0.08997274125
0.09102845458
0.09073768807
0.09085039141
0.0874113711
0.08776610897
0.09017352233
0.08952371679
0.08926983986
0.0889061932
0.08952649284
0.09076799791
0.09215640016
0.08912427254
0.08716605593
0.08961613592
0.09130268951
0.09128187329
0.08790162138
0.08915774108
0.08647174325
0.08540306951
0.08800850308
0.08775034837
0.08844307813
0.08597843204
0.08769072019
0.0843751075
0.07997369538
0.0820785912
0.0813752858
0.08631181899
0.0877320123
0.08589266783
0.0862026795
0.08688474914
0.0839306313
0.08190813164
0.08142688721
0.08044487728
0.08183744319
0.0834769943
0.08237686387
0.0834389161
0.08217677542
0.08162017055
0.08163889287
0.08477527126
0.08551737487
0.08508568443
0.07920828773
0.07998527615
0.07958640999
0.07932811038
0.07805992496
0.07776110338
0.07816797012
0.07977347162
0.07893348477
0.07810760766
0.08336179795
0.08243969882
0.08063561359
0.07845107183
0.07961959968
0.07781141379
0.07773686884
0.08205176301
0.07790285635
0.07285349802
0.07533343234
0.08190924829
0.0814497436
0.0747143943
0.07119830912
0.07052800953
0.07159300455
0.07266799323
0.07491801329
0.07339890869
0.07025002095
0.07465471662
0.07656392459
0.07391202955
0.07122715483
0.06999702788
0.06997300391
0.0720064575
0.07230007924
0.0720182424
0.0729235226
0.0755935729
0.07477588886
0.07070571616
0.06882275185
0.06708629861
0.07189165471
0.07500745539
0.07230376272
0.06859279356
0.06861080337
0.06953960213
0.0722792735
0.06875750442
0.06871920929
0.06718457055
0.06712864023
0.06713863913
0.06548141685
0.06361392089
0.06112019343
0.06212793883
0.06255960398
0.06189451265
0.0616932858
0.06463558735
0.06429455596
0.06358754314
0.0598427313
0.05835709562
0.05823614839
0.05877230807
0.06030957998
0.06150797915
0.06156745461
0.06070614445
0.06031554194
0.05887419981
0.05790346957
0.05511689877
0.05454637964
0.05604405251
0.0591730947
0.0614044369
0.06192214288
0.0576301416
0.0555490907
0.05581644926
0.05806260663
0.05683395503
0.05622826935
0.0559943529
0.05738374585
0.05846787826
0.05740366065
0.0585331412
0.05693748919
0.05810886013
0.05727088865
0.05675404386
0.05510065898
0.05617086193
0.05588792337
0.0574380377
0.05911875362
0.06020375083
0.05869717034
0.05921318401
0.05778074667
0.05679771289
0.05584629762
0.05505985359
0.05595162157
0.05737233143
0.05791389092
0.05735509581
0.05755048533
0.05668233461
0.05589883626
0.05733935235
0.05622605434
0.05575530673
0.05499091731
0.0554843489
0.05770428204
0.05667367873
0.05611213513
0.05524016802
0.05666922283
0.05716826109
0.05841143273
0.05548472809
0.05672374227
0.05583988281
0.05666334558
0.05653233415
0.05486833016
0.05427830771
0.05376663804
0.05435117366
0.05605649823
0.05568418288
0.05440207611
0.05457786352
0.05349669455
0.05352548519
0.05449863867
0.05471505069
0.05492331725
0.05500104066
0.05415394373
0.05303059787
0.05410954034
0.05307933567
0.05386104008
0.05473512483
0.05560919682
0.05434046671
0.054970203
0.05509595895
0.0537580336
0.0550848605
0.05613386291
0.05612604661
0.0550132849
0.05435392175
0.05522991584
0.05492943426
0.05491100172
0.05803776099
0.05661308951
0.05697586901
0.05665447024
0.05690693395
0.05607115655
0.05647917089
0.0555740836
0.05578241573
0.05565301442
0.05462058471
0.05362491982
0.0554719635
0.05623749309
0.05460029446
0.05521248925
0.05669170772
0.0569866044
0.05672919976
0.05598707634
0.05497967026
0.05486621781
0.05509852871
0.05465160456
0.05606551641
0.05561177423
0.05615222768
0.05633265244
0.05589271197
0.05513790438
0.05401290308
0.05436805488
0.05275482305
0.05271531739
0.05376670907
0.05409330036
0.0531529004
0.05353123807
0.05325901261
0.05315042966
0.05297896293
0.05401928919
0.05261069887
0.05441214445
0.05333433455
0.05355231685
0.05320144723
0.05375861551
0.05425245567
0.05400241301
0.05335664836
0.05310370322
0.05370041688
0.0530486745
0.05288076259
0.05374908085
0.0527863796
0.05323651253
0.05475137648
0.05208182889
0.05274864194
0.05280362639
0.05247384137
0.05254684226
0.05385671209
0.05270568531
0.0545190943
0.05418269604
0.05277467856
0.05414824679
0.05314538264
0.05289856517
0.05205194471
0.05251625523
0.05241723612
0.05193693049
0.05316891567
0.05332188765
0.05256892402
0.05132794414
0.05128438894
0.05085427849
0.0499499293
0.0503767354
0.05006616271
0.05025891985
0.05070130767
0.0506606613
0.0509516691
0.0503610574
0.04979193671
0.05011385139
0.04875712144
0.04850124342
0.04910207525
0.04951616582
0.04930273104
0.04975993758
0.04989727889
0.04937950421
0.04844387147
0.04760779641
0.04720146688
0.04710281295
0.04667048126
0.04729997588
0.04743674811
0.0477308222
0.04805525013
0.04729835301
0.04746365559
0.04702814845
0.04590341075
0.04618806766
0.04566628426
0.04559819965
0.04684784155
0.04658240884
0.0465876191
0.04548969666
0.04586737557
0.04472855956
0.04506138745
0.04490940436
0.04560994606
0.04532854581
0.04667353422
0.04636003433
0.04621034576
0.04604108325
0.04576865843
0.0447942238
0.04478830816
0.0448876705
0.04415321151
0.04495949513
0.04505951327
0.04514666452
0.04510988452
0.04440231862
0.04412116496
0.04406009754
0.04337892901
0.04390902096
0.04328683319
0.04358718238
0.0441754927
0.04386167295
0.04438296976
0.04354995084
0.04387082762
0.04350805806
0.04286783855
0.04246115434
0.04254372439
0.04242360032
0.04321173184
0.04282086145
0.04297173562
0.04360909761
0.04276280706
0.04220770429
0.04147156887
0.0413551031
0.04126984525
0.041231596
0.04180793739
0.04185108491
0.04157478508
0.04146304486
0.04135827387
0.04132291203
0.04095721814
0.04098724074
0.04057973935
0.04116436798
0.04095260863
0.04171573739
0.04071025472
0.04116714322
0.04077387541
0.04041230826
0.04019411877
0.04019500066
0.04070894253
0.04033650607
0.03966834526
0.04038739916
0.04024703279
0.04003616631
0.04028432812
0.04000443524
0.04037001721
0.04031983562
0.04007142692
0.04009334979
0.04000723051
0.04010111774
0.04052970469
0.04062382621
0.04095351887
0.03988429571
0.04081498113
0.03984819538
0.04050697301
0.03999020114
0.03957727158
0.03973448472
0.0399936871
0.04011149094
0.03992016593
0.04007128426
0.03946794977
0.03922651096
0.03959247253
0.03939784164
0.03982196025
0.03927515491
0.03952741292
0.03949864564
0.03966172785
0.03927653703
0.03942264187
0.03926511083
0.03926079527
0.03970961173
0.03984082453
0.0402378391
0.04007650523
0.04061212282
0.04053095141
0.03986793623
0.04006734302
0.03985108488
0.03974866266
0.03992185784
0.03991214121
0.04061995369
0.04036690216
0.04045530979
0.04023184694
0.04028766768
0.03980835992
0.03968514552
0.03909013537
0.039600446
0.03985940134
0.0401648902
0.03976651187
0.03938399338
0.03988975755
0.04034393924
0.03994005893
0.0395581485
0.03955595042
0.03973020874
0.03975450838
0.04025861942
0.03997295864
0.04003450898
0.04072272914
0.04063870322
0.03974632374
0.04034482142
0.04056546661
0.04060631353
0.04069603165
0.04065410755
0.04072545046
0.04077296256
0.04052506694
0.04082531089
0.04038810274
0.04050939493
0.04039253133
0.04027195023
0.0405729531
0.04051167053
0.04100799708
0.0411869236
0.04103552917
0.04164855833
0.04125569224
0.04131420948
0.04048316839
0.04021674331
0.04074325889
0.04089506722
0.04100093871
0.04164779825
0.04158666336
0.04158965172
0.04130738705
0.04121278252
0.04123606211
0.04144210551
0.04144917046
0.04164081472
0.041834807
0.04190273745
0.04177071793
0.04153054792
0.0415728496
0.04176495077
0.04163171722
0.04152687484
0.04171420834
0.04135185864
0.04198347323
0.04165260334
0.04215898981
0.04194286669
0.04181611983
0.04163846718
0.04168634853
0.04164382954
0.04142748225
0.0412183763
0.04121124449
0.04171459774
0.04236985134
0.04191346166
0.04177117751
0.04107161241
0.04118791977
0.04122989855
0.04079743322
0.04120190428
0.04113622767
0.04120542924
0.04132986092
0.04107194386
0.04129638209
0.04077732496
0.04100110912
0.04075524474
0.04091986802
0.04055255642
0.04070748067
0.041074936
0.0407165845
0.04063233616
0.04069624982
0.04050432944
0.04020176576
0.04065515013
0.04049067252
0.04064193216
0.0405538437
0.04053546471
0.04038962755
0.0402566898
0.03997371494
0.03938421863
0.03963545158
0.03974291823
0.03981088516
0.0399940971
0.03993916281
0.03952477608
0.03941674628
0.03937465074
0.03965064363
0.03922259301
0.0388210378
0.03941431123
0.03960544485
0.03927495032
0.03913299867
0.03908400998
0.03894718761
0.03872343103
0.03891331797
0.0388501309
0.03871722209
0.03859543688
0.03877526903
0.03889549492
0.03894381287
0.03895674207
0.03911170705
0.03938492594
0.03910803877
0.03876627092
0.03856674942
0.03880619464
0.03870146987
0.03895012746
0.03911606388
0.03913141919
0.03882673973
0.03894692683
0.0389042564
0.03838979495
0.03795017001
0.03796787659
0.03830152423
0.03853531372
0.0387223722
0.03853207099
0.0385870693
0.038440492
0.03868975146
0.03845350798
0.03831401945
0.03875033657
0.03833729193
0.03842511067
0.0384605915
0.03869844312
0.03846132549
0.03882710783
0.03873915736
0.0380677032
0.03825382046
0.03799512446
0.03780206796
0.03802903136
0.0378962264
0.03775947981
0.0382231233
0.0379881708
0.0381020226
0.03794552138
0.03729009495
0.03758124627
0.03772648532
0.03784388724
0.03794717962
0.03829847274
0.03825990418
0.0383581613
0.038239485
0.03835468374
0.03779445097
0.03803863232
0.03817410734
0.03810658884];

t = 1:1:length(rate_ana_1);

figure(1)
plot(t,rate_ana_1*100,'r');
xlabel('tidssteg');
ylabel('Andel anaeroba i procent')
title('Syrenivå 1, 0.2')

hold on
plot(t,rate_ana_02*100,'b');