#include "node.h"

double haversineNode(const Node &a, const Node &b)
{
    const double R = 6371e3;
    double lat1 = a.lat * M_PI / 180;
    double lat2 = b.lat * M_PI / 180;
    double dlat = (b.lat - a.lat) * M_PI / 180;
    double dlon = (b.lon - a.lon) * M_PI / 180;

    double h = sin(dlat / 2) * sin(dlat / 2) +
               cos(lat1) * cos(lat2) * sin(dlon / 2) * sin(dlon / 2);

    return R * 2 * atan2(sqrt(h), sqrt(1 - h));
}

