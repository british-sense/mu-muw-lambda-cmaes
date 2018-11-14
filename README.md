# (*μ*/*μ*<sub>W</sub>, *λ*)-CMA-ES Algorithm


C++における(*μ*/*μ*<sub>W</sub>, *λ*)-CMA-ESの実装です.
実装したCMA-ESは以下の論文のものです.


The CMA Evolution Strategy: A Tutorial
https://hal.inria.fr/hal-01297037/file/tutorial.pdf


以下の環境で実験をしました.

|parameter|value|
|:-:|:-:|
|λ|12|
|μ|6|
|σ|1|
|n|20|
|generation|1000000|
|benchmark|rosenbrock|

実験結果

![rosenbrock](https://github.com/ko-cha/CMA-ES/blob/master/img/frosen.png "rosenbrock")

20次元でのrosenbockにおいて最適解を得ることができました.
