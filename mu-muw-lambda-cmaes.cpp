#include <random>
#include <climits>
#include <algorithm>
#include "algebra.hpp"

using namespace std;
using array = vector<double>;
using matrix = vector<vector<double> >;

random_device seed;
mt19937 mt(0); 
normal_distribution<> dist(0., 1.);

static const int n = 20;
static const double upper_limit = 5;
static const double lower_limit = -5;

class individual{
    
    public:
    array x;
    array z;
    double fitness;

    public:
    individual() : x(zeros(n)), z(zeros(n)), fitness(0.){
    }
    void evaluate(){
        // rastrigin function
        // fitness = 10. * n;
        // for(int i = 0; i < n; i++) fitness += x[i] * x[i] - 10. * cos(2. * M_PI * x[i]);

        // rosenbrock function
        fitness = 0.;
        for(int i = 0; i < n - 1; i++) fitness += 100 * (x[i + 1] - (x[i] * x[i])) * (x[i + 1] - (x[i] * x[i])) + (1 - x[i]) * (1 - x[i]);

        // sphere function
        // fitness = 0.;
        // for(int i = 0; i < n; i++) fitness += x[i] * x[i];

    }
    bool operator<(const individual& rhs)const{
        return fitness < rhs.fitness;
    }
    bool operator>(const individual& rhs)const{
        return fitness > rhs.fitness;
    }
};

matrix extract_col_z(const vector<individual> &indv, const int &start, const int & end){
    matrix mat;
    for(int i = start - 1; i < end; i++){
        mat.push_back(indv[i].z);
    }
    return trans(mat);
}

int main(){

    array x_mean_weight(zeros(n));
    uniform_real_distribution<> range(lower_limit, upper_limit);
    for(int i = 0; i < n; i++) x_mean_weight[i] = range(mt);

    array z_mean_weight(n);
    double sigma = 1.;

    double lambda = 4. + floor(3. * log(n));
    double mu = lambda / 2.;
    array weights = log(mu + (1. / 2.)) - log(seq(1, mu));
    mu = floor(mu);
    weights /= sum(weights);
    double mueff = (sum(weights) * sum(weights)) / sum(weights, 2);

    double c_c = (4. + (mueff / n)) / (n + 4. + (2. * mueff / n));
    double c_s = (mueff + 2.) / (n + mueff + 5.);
    double c_1 = 2. / ((n + 1.3) * (n + 1.3) + mueff);
    double c_mu = min(1. -c_1, 2. * (mueff - 2. + (1. / mueff)) / ((n + 2.) * (n + 2.) + mueff));
    double damp = 1. + (2. * max(0., sqrt((mueff - 1.) / (n + 1.)) - 1.)) + c_s;

    array p_c(zeros(n));
    array p_s(zeros(n));
    matrix B(identity(n));
    matrix D(identity(n));
    matrix C = (B * D) * trans(B * D);
    const double chi_n = sqrt(n) * (1. - (1. / (4. * n)) + (1. / (21. * n * n)));

    double count_evaluate = 0.;
    double stop_evaluate = 1e3 * n * n;
    double stop_fitness = 1e-10;
    double eigen_evaluate = 0.;

    while(count_evaluate < stop_evaluate){

        // 個体の生成
        vector<individual> population(lambda);
        for(int k = 0; k < lambda; k++){
            for(int i = 0; i < n; i++) population[k].z[i] = dist(mt);
            population[k].x = x_mean_weight + sigma * (B * D * population[k].z); 
            population[k].evaluate();
            count_evaluate++;
        }

        // 評価値でソート
        sort(population.begin(), population.end());

        // 平均ベクトルの更新
        x_mean_weight = zeros(x_mean_weight.size());
        z_mean_weight = zeros(z_mean_weight.size());
        for(int k = 0; k < mu; k++){
            x_mean_weight += population[k].x * weights[k];
            z_mean_weight += population[k].z * weights[k];
        }

        // パスの更新
        p_s = (1. - c_s) * p_s + sqrt(c_s * (2. - c_s) * mueff) * (B * z_mean_weight);
        double h_sig = 0.;
        if((norm(p_s) / sqrt(1. - pow((1. - c_s), (2. * count_evaluate / lambda))) / chi_n) < (1.4 + (2. / (n + 1.)))) h_sig = 1.;
        p_c = (1. - c_c) * p_c + h_sig * sqrt(c_c * (2. - c_c) * mueff) * (B * D * z_mean_weight);
        
        // 共分散行列の更新
        C = (1. - c_1 - c_mu) * C 
            + c_1 * (p_c * p_c + (1. - h_sig) * c_c * (2. - c_c) * C)
            + c_mu * (B * D * extract_col_z(population, 1, mu))
            * diagonalization_matrix(weights) * trans(B * D * extract_col_z(population, 1, mu));
        
        // 分散の更新
        sigma *= exp((c_s / damp) * (norm(p_s) / chi_n - 1.));
        
        if((count_evaluate - eigen_evaluate) > (lambda / (c_1 + c_mu) / n / 10.)){
            eigen_evaluate = count_evaluate;
            C = upper_triangle(C) + trans(upper_triangle(C, 1));
            jacobi(C, B, D);
            std::vector<double> dv = sqrt(diagonalization_component(D));
            D = diagonalization_matrix(dv);
        }

        // ステップサイズの調整
        if(population.front().fitness == population[(int)ceil(0.7 * lambda) - 1].fitness) sigma *= exp(0.2 + (c_s / damp));

        // 最良個体の出力
        std::cout << count_evaluate << ',' << population.front().fitness << endl;

        if(population.front().fitness <= stop_fitness) break;
    }
    return 0;
}