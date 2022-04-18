/******************************************************************************
File      : KMEANS_1601.cpp
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
#include <omp.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include "KMeans_1601.h"

using namespace std;

/******************************************************************************
KMeans_1601_Point::arrayToVec()
******************************************************************************/
vector<double> KMeans_1601_Point::arrayToVec(double * vals, int n_vals)
{
    vector<double> values;

    for (int i = 0; i < n_vals; i++)
    {
        values.push_back(vals[i]);
    }

    return values;
}/* end arrayToVec() */

/******************************************************************************
KMeans_1601_Point::lineToVec()
******************************************************************************/
vector<double> KMeans_1601_Point::lineToVec(string &line)
{
    vector<double> values;
    string tmp = "";

    for (int i = 0; i < (int)line.length(); i++)
    {
        if ((48 <= int(line[i]) && int(line[i])  <= 57) || line[i] == '.' || line[i] == '+' || line[i] == '-' || line[i] == 'e')
        {
            tmp += line[i];
        }
        else if (tmp.length() > 0)
        {

            values.push_back(stod(tmp));
            tmp = "";
        }
    }
    if (tmp.length() > 0)
    {
        values.push_back(stod(tmp));
        tmp = "";
    }

    return values;
}/* end lineToVec() */

/******************************************************************************
KMeans_1601_Point()
******************************************************************************/
KMeans_1601_Point::KMeans_1601_Point(int id, string line)
{
    pointId = id;
    values = lineToVec(line);
    dimensions = values.size();
    clusterId = 0; // Initially not assigned to any cluster
}

/******************************************************************************
KMeans_1601_Point()
******************************************************************************/
KMeans_1601_Point::KMeans_1601_Point(int id, double * vals, int n_vals)
{
    pointId = id;
    values = arrayToVec(vals, n_vals);
    dimensions = values.size();
    clusterId = 0; // Initially not assigned to any cluster
}

/******************************************************************************
KMeans_1601_Point::getID()
******************************************************************************/
int KMeans_1601_Point::getID(void) 
{ 
    return pointId; 
}

/******************************************************************************
KMeans_1601_Point::getCluster()
******************************************************************************/
int KMeans_1601_Point::getCluster(void)
{ 
    return clusterId; 
}

/******************************************************************************
KMeans_1601_Point::getDimensions()
******************************************************************************/
int KMeans_1601_Point::getDimensions(void) 
{ 
    return dimensions; 
}

/******************************************************************************
KMeans_1601_Point::setCluster()
******************************************************************************/
void KMeans_1601_Point::setCluster(int val) 
{ 
    clusterId = val; 
}

/******************************************************************************
KMeans_1601_Point::getVal()
******************************************************************************/
double KMeans_1601_Point::getVal(int pos) 
{ 
    return values[pos]; 
}

/******************************************************************************
KMeans_1601_Cluster()
******************************************************************************/
KMeans_1601_Cluster::KMeans_1601_Cluster(int clusterId, KMeans_1601_Point centroid)
{
    this->clusterId = clusterId;
    for (int i = 0; i < centroid.getDimensions(); i++)
    {
        this->centroid.push_back(centroid.getVal(i));
    }
    this->addPoint(centroid);
}

/******************************************************************************
KMeans_1601_Cluster::AddPoint()
******************************************************************************/
void KMeans_1601_Cluster::addPoint(KMeans_1601_Point p)
{
    p.setCluster(this->clusterId);
    points.push_back(p);
}

/******************************************************************************
KMeans_1601_Cluster::removePoint()
******************************************************************************/
bool KMeans_1601_Cluster::removePoint(int pointId)
{
    int size = points.size();

    for (int i = 0; i < size; i++)
    {
        if (points[i].getID() == pointId)
        {
            points.erase(points.begin() + i);
            return true;
        }
    }
    return false;
}

/******************************************************************************
KMeans_1601_Cluster::removeAllPoints()
******************************************************************************/
void KMeans_1601_Cluster::removeAllPoints(void) 
{ 
    points.clear(); 
}

/******************************************************************************
KMeans_1601_Cluster::getId()
******************************************************************************/
int KMeans_1601_Cluster::getId(void) 
{ 
    return clusterId; 
}

/******************************************************************************
KMeans_1601_Cluster::getPoint()
******************************************************************************/
KMeans_1601_Point KMeans_1601_Cluster::getPoint(int pos) 
{ 
    return points[pos]; 
}

/******************************************************************************
KMeans_1601_Cluster::getSize()
******************************************************************************/
int KMeans_1601_Cluster::getSize(void) 
{ 
    return points.size(); 
}

/******************************************************************************
KMeans_1601_Cluster::getCentroidByPos()
******************************************************************************/
double KMeans_1601_Cluster::getCentroidByPos(int pos) 
{ 
    return centroid[pos]; 
}

/******************************************************************************
KMeans_1601_Cluster::setCentroidByPos()
******************************************************************************/
void KMeans_1601_Cluster::setCentroidByPos(int pos, double val) 
{ 
    this->centroid[pos] = val; 
}

/******************************************************************************
KMeans_1601_Alg::clearClusters()
******************************************************************************/
void KMeans_1601_Alg::clearClusters()
{
    for (int i = 0; i < K; i++)
    {
        clusters[i].removeAllPoints();
    }
}

/******************************************************************************
KMeans_1601_Alg::getNearestClusterId()
******************************************************************************/
int KMeans_1601_Alg::getNearestClusterId(KMeans_1601_Point point)
{
    double sum = 0.0, min_dist;
    int NearestClusterId;
    if(dimensions==1) {
        min_dist = abs(clusters[0].getCentroidByPos(0) - point.getVal(0));
    }	
    else 
    {
        for (int i = 0; i < dimensions; i++)
        {
            sum += pow(clusters[0].getCentroidByPos(i) - point.getVal(i), 2.0);
            // sum += abs(clusters[0].getCentroidByPos(i) - point.getVal(i));
        }
        min_dist = sqrt(sum);
    }
    NearestClusterId = clusters[0].getId();

    for (int i = 1; i < K; i++)
    {
        double dist;
        sum = 0.0;
        
        if(dimensions==1) {
                dist = abs(clusters[i].getCentroidByPos(0) - point.getVal(0));
        } 
        else {
            for (int j = 0; j < dimensions; j++)
            {
                sum += pow(clusters[i].getCentroidByPos(j) - point.getVal(j), 2.0);
                // sum += abs(clusters[i].getCentroidByPos(j) - point.getVal(j));
            }

            dist = sqrt(sum);
            // dist = sum;
        }
        if (dist < min_dist)
        {
            min_dist = dist;
            NearestClusterId = clusters[i].getId();
        }
    }

    return NearestClusterId;
}

/******************************************************************************
KMeans_1601_Alg()
******************************************************************************/
KMeans_1601_Alg::KMeans_1601_Alg(int K, int iterations, FILE * pOutfile)
{
    this->K = K;
    this->iters = iterations;
    this->pOut = pOutfile;
}

/******************************************************************************
KMeans_1601_Alg::run()
******************************************************************************/
void KMeans_1601_Alg::run(vector<KMeans_1601_Point> &all_points)
{
    total_points = all_points.size();
    dimensions = all_points[0].getDimensions();

    // Initializing Clusters
    vector<int> used_pointIds;

    for (int i = 1; i <= K; i++)
    {
        while (true)
        {
            int index = rand() % total_points;

            if (find(used_pointIds.begin(), used_pointIds.end(), index) ==
                used_pointIds.end())
            {
                used_pointIds.push_back(index);
                all_points[index].setCluster(i);
                KMeans_1601_Cluster cluster(i, all_points[index]);
                clusters.push_back(cluster);
                break;
            }
        }
    }
    fprintf(pOut, "Clusters initialized = %d\n", clusters.size());
    fprintf(pOut, "Running K-Means Clustering..\n");

    int iter = 1;
    while (true)
    {
        fprintf(pOut, "Iteration %d of %d \n", iter , iters);
        bool done = true;

        // Add all points to their nearest cluster
        #pragma omp parallel for reduction(&&: done) num_threads(16)
        for (int i = 0; i < total_points; i++)
        {
            int currentClusterId = all_points[i].getCluster();
            int nearestClusterId = getNearestClusterId(all_points[i]);

            if (currentClusterId != nearestClusterId)
            {
                all_points[i].setCluster(nearestClusterId);
                done = false;
            }
        }

        // clear all existing clusters
        clearClusters();

        // reassign points to their new clusters
        for (int i = 0; i < total_points; i++)
        {
            // cluster index is ID-1
            clusters[all_points[i].getCluster() - 1].addPoint(all_points[i]);
        }

        // Recalculating the center of each cluster
        for (int i = 0; i < K; i++)
        {
            int ClusterSize = clusters[i].getSize();

            for (int j = 0; j < dimensions; j++)
            {
                double sum = 0.0;
                if (ClusterSize > 0)
                {
                    #pragma omp parallel for reduction(+: sum) num_threads(16)
                    for (int p = 0; p < ClusterSize; p++)
                    {
                        sum += clusters[i].getPoint(p).getVal(j);
                    }
                    clusters[i].setCentroidByPos(j, sum / ClusterSize);
                }
            }
        }

        if (done || iter >= iters)
        {
            fprintf(pOut, "Clustering completed in iteration %d\n", iter);
            break;
        }
        iter++;
    }

    if(pOut != NULL)
    {
        // write cluster assignments to file
        fprintf(pOut, "\n-----------------------------------------------------\n");
        fprintf(pOut, "KMeans cluster assignments\n\n");

        // header
        for (int j = 0; j < dimensions; j++)
        {
            fprintf(pOut,"Param_%d,", j);
        }
        fprintf(pOut,"ClusterID\n");

        // data
        for (int i = 0; i < total_points; i++)
        {
            for (int j = 0; j < dimensions; j++)
            {
                fprintf(pOut,"%f,", all_points[i].getVal(j));
            }
            fprintf(pOut,"%f\n", all_points[i].getCluster());
        }

        // Write cluster centers to file
        fprintf(pOut, "\n-----------------------------------------------------\n");
        fprintf(pOut, "KMeans cluster centers\n\n");

        for (int i = 0; i < K; i++)
        {
            fprintf(pOut, "Cluster # %d centroid : ( ", clusters[i].getId());
            for (int j = 0; j < dimensions; j++)
            {
                fprintf(pOut, "%f ", clusters[i].getCentroidByPos(j));
            }
            fprintf(pOut, ")\n");
        }
    }
}

/******************************************************************************
KMeans_1601_main()

This is the interface function for OSTRICH to use the KMeans_1601 implementation.
******************************************************************************/
int KMeans_1601_main(double ** coords, int n_coords, int n_dims, int k, FILE * pOut)
{
    // Need 3 arguments to run, else exit
    if ((coords == NULL) || (n_coords <= 0) || (n_dims <= 0) || (k <= 0) || (pOut == NULL))
    {
        fprintf(pOut, "KMeans_1601_main() : Error - one or more invlid arguments\n");
        return 1;
    }

    // Fetching number of clusters
    int K = k;

    // fetch points
    int pointId = 1;
    vector<KMeans_1601_Point> all_points;

    for(int i = 0; i < n_coords; i++)
    {
        KMeans_1601_Point point(pointId, coords[i], n_dims);
        all_points.push_back(point);
        pointId++;
    }
    
    fprintf(pOut, "\nData fetched successfully!\n");

    // Return if number of clusters > number of points
    if ((int)all_points.size() < K)
    {
        fprintf(pOut, "Error: Number of clusters greater than number of points.\n");
        return 1;
    }

    // Running K-Means Clustering
    int iters = 100;

    KMeans_1601_Alg kmeans(K, iters, pOut);
    kmeans.run(all_points);

    return 0;
}
