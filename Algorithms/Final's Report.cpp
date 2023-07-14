#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace std;

int M_PI = 3.14;
class CannyEdgeDetector {
public:
    CannyEdgeDetector(const vector<vector<double>>& image4, double sigma = 1.0, double tLow = 0.1, double tHigh = 0.3) :
        image(image4),
        width(image4.size()),
        height(image4[0].size()),
        gradient(width, vector<double>(height)),
        orientation(width, vector<double>(height)),
        nonMaxSuppressed(width, vector<double>(height)),
        visited(width, vector<bool>(height)),
        tLow(tLow),
        tHigh(tHigh) {
        gaussianFilter(sigma);
        computeGradient();
        nonMaximaSuppression();
        doubleThresholding();
    }

    const vector<vector<double>>& getEdges() const {
        return nonMaxSuppressed;
    }

private:
    void gaussianFilter(double sigma) {
        int size = 2 * ceil(2 * sigma) + 1;
        vector<vector<double>> kernel(size, vector<double>(size));
        double sum = 0.0;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                kernel[i][j] = exp(-(pow(i - size / 2, 2) + pow(j - size / 2, 2)) / (2 * pow(sigma, 2)));
                sum += kernel[i][j];
            }
        }
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                kernel[i][j] /= sum;
            }
        }
        for (int i = size / 2; i < width - size / 2; i++) {
            for (int j = size / 2; j < height - size / 2; j++) {
                double val = 0.0;
                for (int k = -size / 2; k <= size / 2; k++) {
                    for (int l = -size / 2; l <= size / 2; l++) {
                        val += image[i + k][j + l] * kernel[size / 2 + k][size / 2 + l];
                    }
                }
                image[i][j] = val;
            }
        }
    }

    void computeGradient() {
        for (int i = 1; i < width - 1; i++) {
            for (int j = 1; j < height - 1; j++) {
                double dx = image[i + 1][j] - image[i - 1][j];
                double dy = image[i][j + 1] - image[i][j - 1];
                gradient[i][j] = sqrt(pow(dx, 2) + pow(dy, 2));
                orientation[i][j] = atan2(dy, dx);
            }
        }
    }

    void nonMaximaSuppression() {
        for (int i = 1; i < width - 1; i++) {
            for (int j = 1; j < height - 1; j++) {
                double theta = orientation[i][j];
                if (theta < 0) theta += M_PI;
                if ((theta >= M_PI * 7 / 8 && theta <= M_PI) || (theta >= 0 && theta <= M_PI / 8)) { // E-W
                    if (gradient[i][j] >= gradient[i][j - 1] && gradient[i][j] >= gradient[i][j + 1]) nonMaxSuppressed[i][j] = gradient[i][j];
                }
                else if (theta >= M_PI / 8 && theta <= M_PI * 3 / 8) { // NE-SW
                    if (gradient[i][j] >= gradient[i - 1][j + 1] && gradient[i][j] >= gradient[i + 1][j - 1]) nonMaxSuppressed[i][j] = gradient[i][j];
                }
                else if (theta >= M_PI * 3 / 8 && theta <= M_PI * 5 / 8) { // N-S
                    if (gradient[i][j] >= gradient[i - 1][j] && gradient[i][j] >= gradient[i + 1][j]) nonMaxSuppressed[i][j] = gradient[i][j];
                }
                else { // NW-SE
                    if (gradient[i][j] >= gradient[i - 1][j - 1] && gradient[i][j] >= gradient[i + 1][j + 1]) nonMaxSuppressed[i][j] = gradient[i][j];
                }
            }
        }
    }

    void doubleThresholding() {
        for (int i = 1; i < width - 1; i++) {
            for (int j = 1; j < height - 1; j++) {
                if (nonMaxSuppressed[i][j] >= tHigh) {
                    nonMaxSuppressed[i][j] = 1.0;
                    dfs(i, j);
                }
            }
        }
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                if (nonMaxSuppressed[i][j] != 1.0) nonMaxSuppressed[i][j] = 0.0;
            }
        }
    }

    void dfs(int i, int j) {
        visited[i][j] = true;
        for (int k = -1; k <= 1; k++) {
            for (int l = -1; l <= 1; l++) {
                if (k == 0 && l == 0) continue;
                int ni = i + k;
                int nj = j + l;
                if (ni >= 0 && ni < width && nj >= 0 && nj < height && !visited[ni][nj]) {
                    if (nonMaxSuppressed[ni][nj] >= tLow) {
                        nonMaxSuppressed[ni][nj] = 1.0;
                        dfs(ni, nj);
                    }
                    else {
                        nonMaxSuppressed[ni][nj] = 0.0;
                    }
                }
            }
        }
    }

    vector<vector<double>> image;
    int width, height;
    vector<vector<double>> gradient;
    vector<vector<double>> orientation;
    vector<vector<double>> nonMaxSuppressed;
    vector<vector<bool>> visited;
    double tLow, tHigh;
};
int main() {
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    vector<vector<double>> image(5, vector<double>(5));
    image[2][2] = 5.0;

    CannyEdgeDetector detector(image);
    const auto& edges = detector.getEdges();
    for (const auto& row : edges) {
        for (double val : row) {
            cout << val << ' ';
        }
        cout << '\n';
    }
    return 0;
}
