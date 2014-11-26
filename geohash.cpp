//
//  geohash.cpp
//
//  Created by Xu, Chen on 14/10/31.
//  Copyright (c) 2014 0d0a.com. All rights reserved.
//

#include <cmath>
#include <algorithm>
#include "geohash.hpp"

////////////////////////////////////////////////////////////////////////////////
// geolocation
////////////////////////////////////////////////////////////////////////////////

/// Global average radii
constexpr long double EARTH_RADIUS=6371.009;
constexpr long double LONGITUDE_CIRCLE=EARTH_RADIUS*2*M_PI;

inline long double radians(long double d) {
    return d * M_PI / 180;
}

inline long double degrees(long double d) {
    return d * 180 / M_PI;
}

inline long double to180(long double d) {
    return d-(std::llround(d)/360*360);
}

inline long double latitude_circle(long double lat) {
    return EARTH_RADIUS*std::cos(radians(lat)) * 2 * M_PI;
}

inline long double longitude_span(long double lat, long double lon1, long double lon2) {
    return radians(std::abs(to180(lon1-lon2)))*latitude_circle(lat)/2/M_PI;
}

inline long double latitude_span(long double lat1, long double lat2) {
    return radians(std::abs(to180(lat1-lat2)))*LONGITUDE_CIRCLE/2/M_PI;
}

double distance(const geolocation &l, const geolocation &r) {
    long double lat1 = radians(l.latitude);
    long double lat2 = radians(r.latitude);
    long double lat_degree = radians(r.latitude - l.latitude);
    long double lon_degree = radians(r.longitude - l.longitude);
    long double a = sin(lat_degree/2.0) * sin(lat_degree/2.0) + cos(lat1) * cos(lat2) * sin(lon_degree /2.0) * sin(lon_degree /2.0);
    long double c = 2.0 * atan2(sqrt(a), sqrt(1.0-a));
    
    return EARTH_RADIUS * c;
}

////////////////////////////////////////////////////////////////////////////////
// bounding_box
////////////////////////////////////////////////////////////////////////////////

bounding_box::bounding_box(geolocation l, double distance) {
    long double latitude_range=degrees(distance/LONGITUDE_CIRCLE);
    long double max_lat=std::max(std::abs(l.latitude-latitude_range), std::abs(l.latitude+latitude_range));
    long double longitude_range=degrees(distance/latitude_circle(max_lat));
    *this=bounding_box(l.latitude-latitude_range,
                       l.latitude+latitude_range,
                       l.longitude-longitude_range,
                       l.longitude+longitude_range);
}

double bounding_box::min_span() const {
    // Returns the shortest side of the box
    return std::min(latitude_span(min_lat, max_lat),
                    std::min(longitude_span(min_lat, min_lon, max_lon),
                             longitude_span(max_lat, min_lon, max_lon)));
}

////////////////////////////////////////////////////////////////////////////////
// geohash
////////////////////////////////////////////////////////////////////////////////

static const char base32_codes[] = {
    '0', '1', '2', '3', '4', '5', '6', '7',
    '8', '9', 'b', 'c', 'd', 'e', 'f', 'g',
    'h', 'j', 'k', 'm', 'n', 'p', 'q', 'r',
    's', 't', 'u', 'v', 'w', 'x', 'y', 'z',
};

static const int base32_indexes[]={
     0,  1,  2,  3,  4,  5,  6,  7, // 30-37, '0'..'7'
     8,  9, -1, -1, -1, -1, -1, -1, // 38-2F, '8','9'
    -1, -1, 10, 11, 12, 13, 14, 15, // 40-47, 'B'..'G'
    16, -1, 17, 18, -1, 19, 20, -1, // 48-4F, 'H','J','K','M','N'
    21, 22, 23, 24, 25, 26, 27, 28, // 50-57, 'P'..'W'
    29, 30, 31, -1, -1, -1, -1, -1, // 58-5F, 'X'..'Z'
    -1, -1, 10, 11, 12, 13, 14, 15, // 60-67, 'b'..'g'
    16, -1, 17, 18, -1, 19, 20, -1, // 68-6F, 'h','j','k','m','n'
    21, 22, 23, 24, 25, 26, 27, 28, // 70-77, 'p'..'w'
    29, 30, 31,                     // 78-7A, 'x'..'z'
};

size_t binary_hash_precision(geolocation l, double dist) {
    for (size_t i = MAX_BINHASH_LENGTH; i >= 1; i--) {
        bounding_box box=decode(binary_encode(l, i));
        if (box.min_span()>=dist*2) {
            return i;
        }
    }
    return 0;
}

binary_hash::binary_hash(const std::string &bit_string)
: bits(0)
, precision(bit_string.size())
{
    for(auto c : bit_string) {
        if (c=='1') {
            bits = (bits<<1) | 1;
        } else {
            bits = bits<<1;
        }
    }
}

binary_hash::binary_hash(geolocation l, double dist) {
    *this=binary_encode(l, binary_hash_precision(l, dist));
}

binary_hash binary_hash::from_geohash(const std::string &hash) {
    binary_hash output;
    for(auto c : hash) {
        int char_index = base32_indexes[c-48];
        if (char_index<0) {
            throw std::invalid_argument("Invalid geohash");
        }
        for (int bits = 4; bits >= 0; --bits) {
            output.push_back(((char_index >> bits) & 1)!=0);
        }
    }
    return output;
}

std::string binary_hash::to_string() const {
    std::string output(precision, ' ');
    uint64_t b=bits<<(64-precision);
    for (int i=0; i<precision; i++) {
        output[i]=(b & 0x8000000000000000LL) ? '1' : '0';
        b<<=1;
    }
    return output;
}

binary_hash binary_encode(geolocation l, size_t precision) {
    // bbox for the lat/lon + errors/ranges
    bounding_box bbox{ -90, 90, -180, 180 };
    bool is_longitude = true;
    
    binary_hash output;
    
    while(output.size() < precision) {
        if (is_longitude) {
            if(l.longitude > bbox.lon_center()) {
                output.push_back(true);
                bbox.min_lon=bbox.lon_center();
            } else {
                output.push_back(false);
                bbox.max_lon=bbox.lon_center();
            }
        } else {
            if(l.latitude > bbox.lat_center() ) {
                output.push_back(true);
                bbox.min_lat = bbox.lat_center();
            } else {
                output.push_back(false);
                bbox.max_lat = bbox.lat_center();
            }
        }
        is_longitude = !is_longitude;
    }
    return output;
}

std::string encode(geolocation l, size_t precision) {
    // DecodedBBox for the lat/lon + errors
    bounding_box bbox{ -90, 90, -180, 180 };
    bool is_longitude = true;
    int num_bits = 0;
    int hash_index = 0;
    
    // Pre-Allocate the hash string
    std::string output(precision, ' ');
    size_t output_length = 0;
    
    while(output_length < precision) {
        if (is_longitude) {
            if(l.longitude > bbox.lon_center()) {
                hash_index = (hash_index << 1) + 1;
                bbox.min_lon=bbox.lon_center();
            } else {
                hash_index = (hash_index << 1) + 0;
                bbox.max_lon=bbox.lon_center();
            }
        } else {
            if(l.latitude > bbox.lat_center() ) {
                hash_index = (hash_index << 1) + 1;
                bbox.min_lat = bbox.lat_center();
            } else {
                hash_index = (hash_index << 1) + 0;
                bbox.max_lat = bbox.lat_center();
            }
        }
        is_longitude = !is_longitude;
        
        ++num_bits;
        if (5 == num_bits) {
            output[output_length] = base32_codes[hash_index];
            output_length++;
            num_bits = 0;
            hash_index = 0;
        }
    }
    return output;
}

bounding_box decode(const binary_hash &hash) {
    // bbox for the lat/lon + errors/ranges
    bounding_box output{ -90, 90, -180, 180 };
    
    bool is_longitude = true;
    
    for(size_t i=1; i<=hash.size(); i++) {
        bool bit = hash.test(i);
        if (is_longitude) {
            if(bit) {
                output.min_lon = output.lon_center();
            } else {
                output.max_lon = output.lon_center();
            }
        } else {
            if(bit) {
                output.min_lat = output.lat_center();
            } else {
                output.max_lat = output.lat_center();
            }
        }
        is_longitude = !is_longitude;
    }
    return output;
}

bounding_box decode(const std::string &hash) {
    bounding_box output{ -90, 90, -180, 180 };
    
    bool is_longitude = true;
    
    for(auto &c : hash) {
        if (c<'0' || c>'z') {
            throw std::invalid_argument("Invalid geohash");
        }
        int char_index = base32_indexes[c-48];
        if (char_index<0) {
            throw std::invalid_argument("Invalid geohash");
        }
        
        for (int bits = 4; bits >= 0; --bits) {
            int bit = (char_index >> bits) & 1;
            if (is_longitude) {
                if(bit == 1) {
                    output.min_lon = output.lon_center();
                } else {
                    output.max_lon = output.lon_center();
                }
            } else {
                if(bit == 1) {
                    output.min_lat = output.lat_center();
                } else {
                    output.max_lat = output.lat_center();
                }
            }
            is_longitude = !is_longitude;
        }
    }
    return output;
}

size_t hash_precision(geolocation l, double dist) {
    for (size_t i = MAX_GEOHASH_LENGTH; i >= 1; i--) {
        bounding_box box=decode(encode(l, i));
        if (box.min_span()>=dist*2) {
            return i;
        }
    }
    return 0;
}

std::string base_hash(geolocation l, double dist) {
    for (size_t i = MAX_GEOHASH_LENGTH; i >= 1; i--) {
        std::string hash=encode(l, i);
        if (decode(hash).min_span()>=dist*2) {
            return hash;
        }
    }
    return "";
}
