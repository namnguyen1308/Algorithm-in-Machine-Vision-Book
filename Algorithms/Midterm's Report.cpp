#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <iostream>
#include <opencv2/opencv.hpp>
using namespace cv;

class ImageProcessor {
private:
    cv::Mat image;

public:
    ImageProcessor(const std::string& imagePath) {
        // Read the image from the specified path
        image = cv::imread(imagePath, cv::IMREAD_GRAYSCALE);

        if (image.empty()) {
            std::cerr << "Failed to open or find the image" << std::endl;
        }
    }

    void applyOpening() {
        int morphSize = 2;
        cv::Mat element = cv::getStructuringElement(
            cv::MORPH_RECT,
            cv::Size(2 * morphSize + 1, 2 * morphSize + 1),
            cv::Point(morphSize, morphSize)
        );

        cv::Mat output;
        cv::morphologyEx(image, output, cv::MORPH_OPEN, element);

        cv::imshow("Original Image", image);
        cv::imshow("Opening Result", output);
        cv::waitKey(0);
    }

    void applyClosing() {
        int morphSize = 2;
        cv::Mat element = cv::getStructuringElement(
            cv::MORPH_RECT,
            cv::Size(2 * morphSize + 1, 2 * morphSize + 1),
            cv::Point(morphSize, morphSize)
        );

        cv::Mat output;
        cv::morphologyEx(image, output, cv::MORPH_CLOSE, element);

        cv::imshow("Original Image", image);
        cv::imshow("Closing Result", output);
        cv::waitKey(0);
    }
};
class IterativeThresholding {
public:
    void segmentImage(const std::string& imagePath) {
        // Read the image
        cv::Mat image = cv::imread(imagePath, cv::IMREAD_GRAYSCALE);
        if (image.empty()) {
            std::cout << "Failed to read the image." << std::endl;
            return;
        }

        // Apply iterative thresholding
        int threshold = getThreshold(image);

        // Segment the image based on the threshold
        cv::Mat segmented;
        cv::threshold(image, segmented, threshold, 255, cv::THRESH_BINARY);

        // Display the segmented image
        cv::imshow("Segmented Image", segmented);
        cv::waitKey(0);
    }

private:
    int getThreshold(const cv::Mat& image) {
        int threshold = 128;  // Initial threshold value
        int prevThreshold = -1;

        // Iterate until the threshold value converges
        while (threshold != prevThreshold) {
            prevThreshold = threshold;

            // Calculate the mean intensity of pixels below and above the threshold
            int sumBelow = 0;
            int countBelow = 0;
            int sumAbove = 0;
            int countAbove = 0;

            for (int i = 0; i < image.rows; i++) {
                for (int j = 0; j < image.cols; j++) {
                    int pixel = image.at<uchar>(i, j);
                    if (pixel < threshold) {
                        sumBelow += pixel;
                        countBelow++;
                    }
                    else {
                        sumAbove += pixel;
                        countAbove++;
                    }
                }
            }

            // Calculate the new threshold as the average of the means
            int meanBelow = (countBelow > 0) ? (sumBelow / countBelow) : 0;
            int meanAbove = (countAbove > 0) ? (sumAbove / countAbove) : 0;
            threshold = (meanBelow + meanAbove) / 2;
        }

        return threshold;
    }
};
class DoubleThresholding {
public:
    DoubleThresholding() {}

    void applyDoubleThresholding(const std::string& imagePath, int lowThreshold, int highThreshold) {
        // Read the image
        cv::Mat image = cv::imread(imagePath, cv::IMREAD_GRAYSCALE);
        if (image.empty()) {
            std::cerr << "Failed to read the image: " << imagePath << std::endl;
            return;
        }

        // Apply double thresholding
        cv::Mat result;
        cv::threshold(image, result, lowThreshold, 255, cv::THRESH_BINARY);
        cv::threshold(result, result, highThreshold, 255, cv::THRESH_BINARY_INV);

        // Display the original and thresholded images
        cv::imshow("Original Image", image);
        cv::imshow("Double Thresholded Image", result);
        cv::waitKey(0);
    }
};
int main() {
//////////////////////////////////////////////////////////////////////
    std::string imagePath = "C:/Users/NAM NGUYEN/Desktop/download.jpg";
    ImageProcessor processor(imagePath);
        processor.applyOpening();
        ////////////////////////
        processor.applyClosing();
 /////////////////////////////////////////////////////////////////////
        IterativeThresholding thresholding;
        thresholding.segmentImage("C:/Users/NAM NGUYEN/Desktop/download.jpg");
//////////////////////////////////////////////////////////////////////
        DoubleThresholding dt;
        dt.applyDoubleThresholding("C:/Users/NAM NGUYEN/Desktop/download.jpg", 100, 200);
//////////////////////////////////////////////////////////////////////
    return 0;
}
