#include "MASTER.h"

using namespace std;

MASTER::MASTER(int K, double delta, int T, int n_init)
    : K(K), delta(delta), T(T), n_init(n_init), c(log(T / delta)), t(1), n(n_init), tn(1), n_hat(log2(T) + 1), alg_index(0) {
    rho = [this](double x) { return sqrt(this->K * this->c / x) + this->K * this->c / x; };
    rho_hat = [this](double x) { return 6 * this->n_hat * this->c * this->rho(x); };
    procedure1();
}

int MASTER::select_arm() {
    int run_len = pow(2, n);
    for (size_t i = 0; i < policy.size(); ++i) {
        if (policy[i].s <= t - tn && t - tn <= policy[i].e && policy[i].len <= run_len) {
            alg_index = i;
            run_len = policy[i].len;
        }
    }

    auto [chosen_arm, f_tilde, theta_hat] = policy[alg_index].select_arm();
    g_tilde.push_back(f_tilde);
    return chosen_arm;
}

void MASTER::update_state(int x, double y) {
    reward.push_back(y);
    policy[alg_index].update_state(x, y);
    int t1=test1();
    int t2=test2();
    t++;
    //If there exists a scenario where MASTER detects any non-stationarity, pause the entire program.
    //Count the detections manually.
    if (t1 == 0 || t2 == 0) {
        cout << "MASTER HAS DETECTED NON-STATIONARITY AT: " << t << ". INITIATING A PAUSE OF THE PROGRAM"<<endl;
        int pause_input; 
        cin >> pause_input;
        n = n_init;
        restart();
    } else {
        if (t >= tn + pow(2, n)) {
            n += 1;
            restart();
        }
    }
}

void MASTER::restart() {
    reward.clear();
    g_tilde.clear();
    active_list.clear();
    tn = t;
    procedure1();
}

void MASTER::re_init() {
    t = 1;
    n = n_init;
    tn = t;
    policy.clear();
    alg_index = 0;
    g_tilde.clear();
    alge.clear();
    reward.clear();
    procedure1();
}

void MASTER::procedure1() {
    policy.clear();
    for (int tau = 0; tau < pow(2, n); ++tau) {
        for (int i = 0; i <= n; ++i) {
            int m = n - i;
            if (tau % int(pow(2, m)) == 0 && distribution(generator) < rho(pow(2, n)) / rho(pow(2, m))) {
                int algs = tau;
                int alge = tau + pow(2, m) - 1;
                policy.emplace_back(UCB(T, delta, algs, alge));
            }
        }
    }
}

int MASTER::test1() {
    double U = *min_element(g_tilde.begin(), g_tilde.end());
    for (const auto& alg : policy) {
        if (t - tn == alg.e) {
            double R_sum = accumulate(reward.begin() + alg.s, reward.begin() + alg.e + 1, 0.0);
            if (R_sum / alg.len >= U + 9 * rho_hat(alg.len)) {
                return 0;
            }
        }
    }
    return 1;
}

int MASTER::test2() {
    double a_sum = accumulate(g_tilde.begin(), g_tilde.end(), 0.0) - accumulate(reward.begin(), reward.end(), 0.0);
    if (a_sum / (t - tn + 1) >= 3 * rho_hat(t - tn + 1)) {
        return 0;
    }
    return 1;
}

vector<int> MASTER::run(const vector<vector<int>>& Environment) {
    vector<int> ChosenArms;
    for (int i = 0; i < T; ++i) {
        int chosen_arm = select_arm();
        ChosenArms.push_back(chosen_arm);
        update_state(chosen_arm, Environment[chosen_arm][t - 1]);
    }
    return ChosenArms;
}
