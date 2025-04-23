#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include "helpers.h"
#include "../TFTFGraph.h"
#include <algorithm>

Coordinate interpolate(const Coordinate &start, const Coordinate &end, float t)
{
    float lat = start.latitude + t * (end.latitude - start.latitude);
    float lon = start.longitude + t * (end.longitude - start.longitude);
    return {lat, lon};
}

std::vector<Coordinate> densifyPath(const std::vector<Coordinate> &path, float spacingMeters)
{
    std::vector<Coordinate> result;
    for (size_t i = 1; i < path.size(); ++i)
    {
        const auto &start = path[i - 1];
        const auto &end = path[i];
        float dist = haversine(start, end);
        int numPoints = static_cast<int>(dist / spacingMeters);

        if (numPoints == 0)
        {
            result.push_back(start); // optional: skip if you don't want to repeat start
        }
        else
        {
            for (int j = 0; j < numPoints; ++j)
            {
                float t = static_cast<float>(j) / numPoints;
                result.push_back(interpolate(start, end, t));
            }
        }
    }
    result.push_back(path.back());
    return result;
}


float haversine(const Coordinate &a, const Coordinate &b)
{
    auto deg2rad = [](double deg)
    {
        return deg * M_PI / 180.0;
    };

    double lat1 = deg2rad(a.latitude);
    double lon1 = deg2rad(a.longitude);
    double lat2 = deg2rad(b.latitude);
    double lon2 = deg2rad(b.longitude);

    double dlat = lat2 - lat1;
    double dlon = lon2 - lon1;

    double h = std::sin(dlat / 2) * std::sin(dlat / 2) +
               std::cos(lat1) * std::cos(lat2) *
                   std::sin(dlon / 2) * std::sin(dlon / 2);

    double c = 2 * std::atan2(std::sqrt(h), std::sqrt(1 - h));
    return EARTH_RADIUS_METERS * c;
}

int closestCoordinateIndex(const std::vector<Coordinate> &path, const Coordinate &target)
{
    int bestIdx = 0;
    float bestDist = std::numeric_limits<float>::infinity();
    for (size_t i = 0; i < path.size(); ++i)
    {
        float dist = haversine(path[i], target);
        if (dist < bestDist)
        {
            bestDist = dist;
            bestIdx = i;
        }
    }
    return bestIdx;
}

float computeRouteDistance(const std::vector<Coordinate> &path, const Coordinate &start, const Coordinate &end)
{
    if (path.empty())
        return 0.0f;

    int closestStart = -1;
    int closestEnd = -1;
    float minStartDist = std::numeric_limits<float>::infinity();
    float minEndDist = std::numeric_limits<float>::infinity();

    for (int i = 0; i < path.size(); ++i)
    {
        float dStart = haversine(path[i], start);
        float dEnd = haversine(path[i], end);
        if (dStart < minStartDist)
        {
            minStartDist = dStart;
            closestStart = i;
        }
        if (dEnd < minEndDist)
        {
            minEndDist = dEnd;
            closestEnd = i;
        }
    }

    if (closestStart > closestEnd)
        std::swap(closestStart, closestEnd);

    float distance = 0.0f;
    for (int i = closestStart; i < closestEnd; ++i)
    {
        distance += haversine(path[i], path[i + 1]);
    }

    return distance;
}

double getActualSegmentDistance(const Coordinate &start, const Coordinate &end, const std::vector<Coordinate> &path, bool isLoop)
{
    int startIndex = getClosestIndex(path, start);
    int endIndex = getClosestIndex(path, end);

    double distance = 0.0;

    // Prevent illegal wraparound
    if (!isLoop && startIndex > endIndex)
    {
        std::cerr << "Warning: Non-loop route but segment wraps around. Returning large distance.\n";
        return 1e9;
    }

    if (startIndex <= endIndex)
    {
        for (int i = startIndex; i < endIndex; ++i)
        {
            distance += haversine(path[i], path[i + 1]);
        }
    }
    else
    { // loop case
        for (int i = startIndex; i < path.size() - 1; ++i)
        {
            distance += haversine(path[i], path[i + 1]);
        }
        distance += haversine(path.back(), path[0]);
        for (int i = 0; i < endIndex; ++i)
        {
            distance += haversine(path[i], path[i + 1]);
        }
    }

    return distance;
}

float getSubpathDistance(const std::vector<Coordinate> &coords, int i1, int i2, bool isLoop)
{
    float total = 0.0f;
    if (i1 == i2)
        return 0.0f;

    if (!isLoop || i1 < i2)
    {
        for (int i = i1; i < i2; ++i)
        {
            total += haversine(coords[i], coords[i + 1]);
        }
    }
    else
    {
        for (int i = i1; i + 1 < coords.size(); ++i)
        {
            total += haversine(coords[i], coords[i + 1]);
        }
        for (int i = 0; i < i2; ++i)
        {
            total += haversine(coords[i], coords[i + 1]);
        }
    }
    return total;
}

double roundUpToNearest2_5(double amount)
{
    return std::ceil(amount / 2.5) * 2.5;
}

Coordinate projectOntoPath(const Coordinate &point, const std::vector<Coordinate> &path)
{
    double minDistance = std::numeric_limits<double>::max();
    Coordinate closestPoint;

    for (size_t i = 0; i + 1 < path.size(); ++i)
    {
        const Coordinate &A = path[i];
        const Coordinate &B = path[i + 1];

        double dx = B.longitude - A.longitude;
        double dy = B.latitude - A.latitude;

        double lengthSquared = dx * dx + dy * dy;
        if (lengthSquared == 0.0)
            continue;

        double t = ((point.longitude - A.longitude) * dx + (point.latitude - A.latitude) * dy) / lengthSquared;
        t = std::max(0.0, std::min(1.0, t));

        Coordinate projection = {
            A.latitude + t * dy,
            A.longitude + t * dx};

        double d = haversine(point, projection);
        if (d < minDistance)
        {
            minDistance = d;
            closestPoint = projection;
        }
    }

    return closestPoint;
}

int getClosestIndex(const std::vector<Coordinate> &path, const Coordinate &coord)
{
    float minDist = std::numeric_limits<float>::infinity();
    int index = -1;

    for (size_t i = 0; i < path.size(); ++i)
    {
        float dist = haversine(coord, path[i]);
        if (dist < minDist)
        {
            minDist = dist;
            index = static_cast<int>(i);
        }
    }

    return index;
}

std::vector<Coordinate> getFullSegmentPath(const std::vector<Coordinate> &path, int startIdx, int endIdx, bool isLoop)
{
    std::vector<Coordinate> segmentPath;

    if (isLoop && startIdx > endIdx)
    {
        // Wrap-around: insert from startIdx to end, then from 0 to endIdx
        segmentPath.insert(segmentPath.end(), path.begin() + startIdx, path.end());
        segmentPath.insert(segmentPath.end(), path.begin(), path.begin() + endIdx + 1);
    }
    else
    {
        // Normal slice
        segmentPath.insert(segmentPath.end(), path.begin() + startIdx, path.begin() + endIdx + 1);
    }

    return segmentPath;
}

std::vector<int> getCandidatePathIndexes(
    const std::vector<Coordinate> &path,
    const Coordinate &target,
    float radiusMeters = 30.0f,
    int minIndexSeparation = 5)
{
    std::vector<int> candidates;

    for (size_t i = 0; i < path.size(); ++i)
    {
        float dist = haversine(path[i], target);
        if (dist <= radiusMeters)
        {
            // Skip if index is too close to any existing candidate
            bool tooClose = false;
            for (int c : candidates)
            {
                if (std::abs(static_cast<int>(i) - c) < minIndexSeparation)
                {
                    tooClose = true;
                    break;
                }
            }
            if (!tooClose)
                candidates.push_back(static_cast<int>(i));
        }
    }

    return candidates;
}

std::vector<Coordinate> getShortestSegmentPath(
    const std::vector<Coordinate> &path,
    const Coordinate &startCoord,
    const Coordinate &endCoord,
    bool isLoop)
{
    const float radiusMeters = 50.0f;
    const int minIndexSeparation = 5;

    std::vector<int> startCandidates = {getClosestIndex(path, startCoord)};
    std::vector<int> endCandidates = getCandidatePathIndexes(path, endCoord, radiusMeters, minIndexSeparation);

    // Fallback: add closest index if none found
    if (startCandidates.empty())
        startCandidates.push_back(getClosestIndex(path, startCoord));
    if (endCandidates.empty())
        endCandidates.push_back(getClosestIndex(path, endCoord));

    std::vector<Coordinate> bestPath;
    float bestDist = std::numeric_limits<float>::infinity();

    for (int startIdx : startCandidates)
    {
        for (int endIdx : endCandidates)
        {
            std::vector<Coordinate> forwardPath;
            std::vector<Coordinate> loopPath;

            if (startIdx <= endIdx)
            {
                forwardPath.insert(forwardPath.end(), path.begin() + startIdx, path.begin() + endIdx + 1);
            }
            else
            {
                forwardPath.insert(forwardPath.end(), path.begin() + startIdx, path.end());
                forwardPath.insert(forwardPath.end(), path.begin(), path.begin() + endIdx + 1);
            }

            float forwardDist = 0.0f;
            for (size_t i = 1; i < forwardPath.size(); ++i)
                forwardDist += haversine(forwardPath[i - 1], forwardPath[i]);

            float loopDist = std::numeric_limits<float>::infinity();
            if (isLoop)
            {
                if (endIdx <= startIdx)
                {
                    loopPath.insert(loopPath.end(), path.begin() + endIdx, path.begin() + startIdx + 1);
                }
                else
                {
                    loopPath.insert(loopPath.end(), path.begin() + endIdx, path.end());
                    loopPath.insert(loopPath.end(), path.begin(), path.begin() + startIdx + 1);
                }
                std::reverse(loopPath.begin(), loopPath.end());

                loopDist = 0.0f;
                for (size_t i = 1; i < loopPath.size(); ++i)
                    loopDist += haversine(loopPath[i - 1], loopPath[i]);
            }

            // Pick the better path from forward vs loop for this candidate pair
            if (forwardDist <= loopDist && forwardDist < bestDist)
            {
                bestDist = forwardDist;
                bestPath = forwardPath;
            }
            else if (isLoop && loopDist < forwardDist && loopDist < bestDist)
            {
                bestDist = loopDist;
                bestPath = loopPath;
            }
        }
    }

    return bestPath;
}
