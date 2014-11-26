//
//  geohash.hpp
//
//  Created by Xu, Chen on 14/10/31.
//  Copyright (c) 2014 0d0a.com. All rights reserved.
//

#ifndef geohash_hpp_included
#define geohash_hpp_included

#include <algorithm>
#include <string>

constexpr size_t MAX_GEOHASH_LENGTH=12;
constexpr size_t MAX_BINHASH_LENGTH=64;

/// WGS84 point
struct geolocation {
    double latitude;
    double longitude;
};

inline bool operator==(const geolocation &l, const geolocation &r)
{ return (l.latitude==r.latitude) && (l.longitude==r.longitude); }

inline bool operator!=(const geolocation &l, const geolocation &r)
{ return (l.latitude!=r.latitude) || (l.longitude!=r.longitude); }

/// Distance of 2 geolocations
double distance(const geolocation &l, const geolocation &r);

inline double operator-(const geolocation &l, const geolocation &r)
{ return distance(l, r); }

/// Latitude/longitude bounding box
struct bounding_box {
    bounding_box()=default;
    
    bounding_box(double lat1, double lat2, double lon1, double lon2)
    : min_lat(std::min(lat1, lat2))
    , max_lat(std::max(lat1, lat2))
    , min_lon(std::min(lon1, lon2))
    , max_lon(std::max(lon1, lon2))
    {}
    
    bounding_box(geolocation l1, geolocation l2)
    : min_lat(std::min(l1.latitude, l2.latitude))
    , max_lat(std::max(l1.latitude, l2.latitude))
    , min_lon(std::min(l1.longitude, l2.longitude))
    , max_lon(std::max(l1.longitude, l2.longitude))
    {}
    
    /// Create a bounding box that contains the circle
    bounding_box(geolocation l, double distance);
    
    /// Test if a point in the box
    bool contains(geolocation l) const {
        return (l.latitude>=min_lat) && (l.latitude<=max_lat) && (l.longitude>=min_lon) && (l.longitude<=max_lon);
    }
    
    /// Center
    geolocation center() const { return geolocation{lat_center(), lon_center()};}
    /// Corners
    geolocation bottom_left() const { return geolocation{min_lat, min_lon}; }
    geolocation bottom_right() const { return geolocation{min_lat, max_lon}; }
    geolocation top_left() const { return geolocation{max_lat, min_lon}; }
    geolocation top_right() const { return geolocation{max_lat, max_lon}; }
    
    /// Edges
    double lat_range() const { return max_lat-min_lat; }
    double lon_range() const { return max_lon-min_lon; }
    double lat_err() const { return (max_lat-min_lat)/2; }
    double lon_err() const { return (max_lon-min_lon)/2; }
    double lat_center() const { return (min_lat+max_lat)/2; }
    double lon_center() const { return (min_lon+max_lon)/2; }

    /// Merge with other box
    bounding_box &merge(const bounding_box &b) {
        min_lat=std::min(min_lat, b.min_lat);
        max_lat=std::max(max_lat, b.max_lat);
        min_lon=std::min(min_lon, b.min_lon);
        max_lon=std::max(max_lon, b.max_lon);
        return *this;
    }
    
    /// Returns the shortest edge of the box
    double min_span() const;
    
    /// Default value is the largest box
    double min_lat=-90.0;
    double max_lat=90.0;
    double min_lon=-180.0;
    double max_lon=180.0;
};

inline bool operator==(const bounding_box &b1, const bounding_box &b2) {
    return (b1.min_lat==b2.min_lat) && (b1.max_lat==b2.max_lat)
    && (b1.min_lon==b2.min_lon) && (b1.max_lon==b2.max_lon);
}

inline bool operator!=(const bounding_box &b1, const bounding_box &b2) {
    return !(b1==b2);
}

/// Merge 2 bounding boxes
inline bounding_box merge(const bounding_box &b1, const bounding_box &b2) {
    return bounding_box{
        std::min(b1.min_lat, b2.min_lat),
        std::max(b1.max_lat, b2.max_lat),
        std::min(b1.min_lon, b2.min_lon),
        std::max(b1.max_lon, b2.max_lon),
    };
}

/// Binary hash code
struct binary_hash {
    binary_hash()=default;
    binary_hash(uint64_t b, size_t p) : bits(b), precision(p) {}
    binary_hash(const std::string &bit_string);
    binary_hash(geolocation l, double dist);

    static binary_hash from_geohash(const std::string &hash);
    
    size_t size() const { return precision; }
    bool empty() const { return size()==0; }
    bool test(size_t n) const { return (bits & (1ull << (precision-n)))!=0; }
    void push_back(bool b) { bits<<=1; bits|=(b?1:0); precision++; }

    operator std::string() const { return to_string(); }
    
    std::string to_string() const;

    uint64_t bits=0;
    size_t precision=0;
};

inline bool operator==(const binary_hash &b1, const binary_hash &b2) {
    return b1.empty() ? b2.empty() : (b1.bits==b2.bits && b1.precision==b2.precision);
}

inline bool operator!=(const binary_hash &b1, const binary_hash &b2) {
    return !(b1==b2);
}

/// Required precision to represent the area
size_t binary_hash_precision(geolocation l, double dist);
/// Binary encode with specific precision
binary_hash binary_encode(geolocation l, size_t bit_count);
/// Decode binary code into a bounding box
bounding_box decode(const binary_hash &bits);
/// Get the neighbor on specific direction of this binary code
inline binary_hash neighbor(const binary_hash &hash,
                            const std::pair<int, int> &direction)
{
    bounding_box b=decode(hash);
    geolocation cp=b.center();
    cp.latitude += direction.first * b.lat_range();
    cp.longitude += direction.second * b.lon_range();
    return binary_encode(cp, hash.size());
}

/// Base32 hash string
std::string encode(geolocation l, size_t precision);
bounding_box decode(const std::string &hash);
inline std::string neighbor(const std::string &hash,
                            const std::pair<int, int> &direction)
{
    bounding_box b=decode(hash);
    geolocation cp=b.center();
    cp.latitude += direction.first * b.lat_range();
    cp.longitude += direction.second * b.lon_range();
    return encode(cp, hash.size());
}

/// Encode with specific ranges
template<typename Container>
void encode_precision_range(geolocation l,
                            std::back_insert_iterator<Container> i,
                            size_t range_largest=1,
                            size_t range_smallest=MAX_GEOHASH_LENGTH)
{
    for (size_t n=range_largest; n<range_smallest+1; n++) {
        *i++=encode(l, n);
    }
}

/// Encode with specific ranges
template<typename Container>
void encode_precision_range(geolocation l,
                            Container &c,
                            size_t range_largest=1,
                            size_t range_smallest=MAX_GEOHASH_LENGTH)
{
    encode_precision_range(l, std::back_inserter(c), range_largest, range_smallest);
}

/// Test if the point in the bounding box represented by this hash code
inline bool hash_contains(const std::string &hash, geolocation l) {
    return decode(hash).contains(l);
}

/// Returns hash precision for given location and distance range
size_t hash_precision(geolocation l, double dist);

/// Returns hash code for the smallest bounding box contains the location and has span longer than dist*2
/// The circle range is contained by this box and its neighbors
/// Use this instead of one big box because one big box may actually contain 32 smaller boxes
std::string base_hash(geolocation l, double dist);

/// Get the minimal geohash and its neighbors for given location and range
/// Returns 9 geohash codes instead of one big box, which contains 32 smaller boxes
template<typename Container>
void hash_codes(geolocation l, double dist, std::back_insert_iterator<Container> i) {
    std::string hash=base_hash(l, dist);
    *i++=neighbor(hash, {-1, -1});
    *i++=neighbor(hash, {-1,  0});
    *i++=neighbor(hash, {-1,  1});
    *i++=neighbor(hash, { 0, -1});
    *i++=neighbor(hash, { 0,  1});
    *i++=neighbor(hash, { 1, -1});
    *i++=neighbor(hash, { 1,  0});
    *i++=neighbor(hash, { 1,  1});
    *i++=std::move(hash);
}

template<typename Container>
void hash_codes(geolocation l, double dist, Container &c) {
    hash_codes(l, dist, std::back_inserter(c));
}

#endif
