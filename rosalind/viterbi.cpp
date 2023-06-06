// ba10c

#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <set>
#include <map>

using namespace std;

typedef long long int ll;

int main() {
    string div = "--------";
    string tmp;

    string seq;
    cin >> seq;
    cin >> tmp;
    map<string, int> v_to_i;
    map<int, string> i_to_v;
    int size_v = 0;
    while (true) {
        cin >> tmp;
        if (tmp == div) break;
        int index = size_v++;
        v_to_i[tmp] = index;
        i_to_v[index] = tmp;
    }


    map<string, int> h_to_i;
    map<int, string> i_to_h;
    int size_h = 0;
    while (true) {
        cin >> tmp;
        if (tmp == div) break;
        int index = size_h++;
        h_to_i[tmp] = index;
        i_to_h[index] = tmp;
    }

    vector<vector<double>> matrix(size_h, vector<double>(size_h));
    vector<int> columns_h(size_h);
    for (int i = 0; i < size_h; ++i) {
        string col;
        cin >> col;
        columns_h[i] = h_to_i[col];
    }
    for (int i = 0; i < size_h; ++i) {
        string row;
        cin >> row;
        int row_i = h_to_i[row];
        for (int j = 0; j < size_h; ++j) {
            double cell;
            cin >> cell;
            matrix[row_i][columns_h[j]] = cell;
        }
    }

    cin >> tmp;

    vector<vector<double>> emission(size_h, vector<double>(size_v));
    vector<int> columns_v(size_v);
    for (int i = 0; i < size_v; ++i) {
        string col;
        cin >> col;
        columns_v[i] = v_to_i[col];
    }
    for (int i = 0; i < size_h; ++i) {
        string row;
        cin >> row;
        int row_i = h_to_i[row];
        for (int j = 0; j < size_v; ++j) {
            double cell;
            cin >> cell;
            emission[row_i][columns_v[j]] = cell;
        }
    }

    vector<int> num_seq;
    for (char c : seq) {
        num_seq.push_back(v_to_i[string(1, c)]);
    }
    vector<vector<double>> probs(num_seq.size(), vector<double>(size_h));
    vector<vector<int>> prev(num_seq.size(), vector<int>(size_h, -1));
    for (int i = 0; i < num_seq.size(); ++i) {
        int cur = num_seq[i];
        if (i == 0) {
            for (int j = 0; j < size_h; ++j) {
                probs[i][j] = emission[j][cur];
            }
        } else {
            for (int j = 0; j < size_h; ++j) {
                double best_prob = 0;
                int best_prev = -1;
                for (int k = 0; k < size_h; ++k) {
                    double prob = emission[j][cur] * matrix[k][j] * probs[i - 1][k];
                    if (prob > best_prob) {
                        best_prob = prob;
                        best_prev = k;
                    }
                }
                probs[i][j] = best_prob;
                prev[i][j] = best_prev;
            }
        }
    }
    double max_finish = 0;
    int max_finish_i = 0;
    int s_i = num_seq.size() - 1;
    for (int i = 0; i < size_h; ++i) {
        double cur = probs[s_i][i];
        if (cur > max_finish) {
            max_finish = cur;
            max_finish_i = i;
        }
    }
    vector<string> rev_ans;
    while (max_finish_i != -1) {
        rev_ans.push_back(i_to_h[max_finish_i]);
        max_finish_i = prev[s_i][max_finish_i];
        --s_i;
    }
    for (int i = rev_ans.size() - 1; i >= 0; --i) {
        cout << rev_ans[i];
    }
    return 0;
}