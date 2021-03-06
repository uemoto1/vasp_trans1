# vasp_trans1

__電気双極子近似に基づく状態間遷移確率計算プログラム__

このプログラムでは、VASP出力のWAVECARファイルに含まれる擬波動関数データをもとに、
ガンマ点上(k=0)の２つのバンド間の遷移確率計算を行う。
また必要に応じてシザーズ補正された計算を実行できる。

## ビルド方法
* 標準的なCコンパイラでは以下のコマンドライン引数でコンパイル可能
 * gcc(Gnu C Compiler) 4.8以降
    
    gcc trans1.c -o trans1 -std=c99 -lm -O3

 * icc(Intel C Compiler)
    
    icc trans1.c -o trans1 -std=c99 -lm -O3

## 使用方法
* WAVECARファイルを含むディレクトリ上でプログラムを実行する

        trans1 ~/a_directory_containing_wavecar_file/WAVECAR

* またScissors Approximationによるスペクトルの補正を行う場合

        trans1 ~/a_directory_containing_wavecar_file/WAVECAR 1.00 0.00
        
        trans1 [INPUTFILE] [SCISSORS DELTA] [SCISSORS_OFFSET]
        
    * [SCISSORS DELTA]
        * シザーズ補正の大きさ
    * [SCISSORS OFFSET]
        * 補正の挿入位置（フェルミ準位を基準にする）


* 結果は標準出力に書きだされる

        # eigen energies on k=0
        EIGEN i=0 e=-36.577988 n=1.000000
        EIGEN i=1 e=-36.390068 n=1.000000
        （中略）
        # transition at k=0
        # i=initial state index
        # j=excited state index
        # e=energy difference between two states
        # s=scissors approximated energy difference
        # n=population difference between two states
        # wx, wy, wr, wl = transition probability (x,y,l,r-polarization)
        TRANS i=798 j=957 e=2.499406 s=3.499406 n=0.510406 wx=2.878050e+03 wy=3.314375e-01 wl=3.754002e+04 wr=7.269572e+04
        TRANS i=799 j=957 e=2.487174 s=3.487174 n=0.510406 wx=2.711653e+05 wy=7.709108e+04 wl=2.616695e+03 wr=6.202896e+05

* EIGEN: 固有状態
    * i=バンド番号 e=エネルギー n=占有率
* TRANS 遷移確率
    * i=始状態 j=終状態 e=エネルギー差 s=シザーズ補正されたエネルギー n=占有率差
    * wx, wy, wl, wr : 直線偏光(X, Y)および円偏光(L, R)による遷移確率

## 更新履歴
* シザーズ補正用ルーチンを追加
* コマンドライン引数による入力ファイルの指定、Usageの表示を追加

