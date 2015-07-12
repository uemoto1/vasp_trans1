# vasp_trans1

__電気双極子近似に基づく状態間遷移確率計算プログラム__

このプログラムでは、VASP出力のWAVECARファイルに含まれる擬波動関数データをもとに、
ガンマ点上(k=0)の２つのバンド間の遷移確率計算を行う。

## ビルド方法
* 標準的なCコンパイラでは以下のコマンドライン引数でコンパイル可能

    gcc trans1.c -o trans1 -std=c99 -lm -O3

## 使用方法
* WAVECARファイルを含むディレクトリ上でプログラムを実行する

        cd ~/directory_contains_vasp_output/
        ~/directory_of_trans1/trans1

* 出力ファイルは以下のようになる

        # eigen energies on k=0
        EIGEN i=0 e=-36.577988 n=1.000000
        EIGEN i=1 e=-36.390068 n=1.000000
        （中略）
        # transition at k=0
        # i=initial state index
        # j=excited state index
        # e=energy difference between two states
        # n=population difference between two states
        # wx, wy, wr, wl = transition probability (x,y,l,r-polarization)
        TRANS i=798 j=957 e=2.499406 n=0.510406 wx=2.878050e+03 wy=3.314375e-01 wl=3.754002e+04 wr=7.269572e+04
        TRANS i=799 j=957 e=2.487174 n=0.510406 wx=2.711653e+05 wy=7.709108e+04 wl=2.616695e+03 wr=6.202896e+05

* EIGEN: 固有状態
    * i=バンド番号 e=エネルギー n=占有率
* TRANS 遷移確率
    * i=始状態 j=終状態 e=エネルギー差 n=占有率差
    * wx, wy, wl, wr : 直線偏光(X, Y)および円偏光(L, R)による遷移確率
