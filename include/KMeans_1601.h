/******************************************************************************
File      : KMEANS_1601.h
Author    : Aditya Sharma and L. Shawn Matott 
Copyright : 2022, Aditya Sharma and L. Shawn Matott

This is a C++ implementation of the simple K-Means clustering algorithm.

It is available on github (https://github.com/aditya1601/kmeans-clustering-cpp) 
and distributed via the MIT license. It was developed by Aditya Sharma and adapted 
for use in OSTRICH by L. Shawn Matott.

K-means clustering is a type of unsupervised learning, which is used when you have 
unlabeled data (i.e., data without defined categories or groups). The goal of this 
algorithm is to find groups in the data, with the number of groups represented by 
the variable K. The algorithm works iteratively to assign each data point to one of 
K groups based on the features that are provided. Data points are clustered based 
on feature similarity.

Version History
04-18-2022   lsm   added copyright information and initial comments.
******************************************************************************/
#ifndef KMEANS_1601_H
#define KMEANS_1601_H

#include <vector>
#include <string>

class KMeans_1601_Point
{
    public:
        KMeans_1601_Point(int id, std::string line);
        KMeans_1601_Point(int id, double * vals, int n_vals);
        int getDimensions(void);
        int getCluster(void);
        int getID(void);
        void setCluster(int val);
        double getVal(int pos);

    private:
        int pointId;
        int clusterId;
        int dimensions;
        std::vector<double> values;
        std::vector<double> lineToVec(std::string &line);
        std::vector<double> arrayToVec(double * vals, int n_vals);
};

class KMeans_1601_Cluster
{
    public:
        KMeans_1601_Cluster(int clusterId, KMeans_1601_Point centroid);
        void addPoint(KMeans_1601_Point p);
        bool removePoint(int pointId);
        void removeAllPoints(void);
        int getId(void);
        KMeans_1601_Point getPoint(int pos);
        int getSize(void);
        double getCentroidByPos(int pos);
        void setCentroidByPos(int pos, double val);

    private:
        int clusterId;
        std::vector<double> centroid;
        std::vector<KMeans_1601_Point> points;
};

class KMeans_1601_Alg
{
    public:
        KMeans_1601_Alg(int K, int iterations, FILE * pOutfile);
        void run(std::vector<KMeans_1601_Point> &all_points, double ** centers);

    private:
        int K, iters, dimensions, total_points;
        std::vector<KMeans_1601_Cluster> clusters;
        FILE * pOut;

        void clearClusters(void);
        int getNearestClusterId(KMeans_1601_Point point);
};

double ** KMeans_1601_main(double ** coords, int n_coords, int n_dims, int k, FILE * pOut);

#endif /* KMEANS_1601_H */
