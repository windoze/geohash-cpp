#include <vector>
#include <assert.h>
#include "geohash.hpp"

void test_geolocation() {
	geolocation l1{12.5, 26.75};
	geolocation l2{22.5, 88.75};
	geolocation l3{12.5, 26.75};
	assert(l1!=l2);
	assert(l1==l3);
	// Distance is about 6625.8km
	assert(((l1-l2)>6625) && ((l1-l2)<6626));
}

// min/max correction
void test_bbox1() {
	bounding_box b1{100, 50, 20, 40};
	bounding_box b2{50, 100, 20, 40};
	assert(b1==b2);
}

// Merge
void test_bbox2() {
	bounding_box b1{12, 34, 56, 78};
	bounding_box b2{34, 56, 78, 90};
	bounding_box b3{12, 56, 56, 90};
	bounding_box b4{13, 35, 57, 79};
	assert(b3==merge(b1, b2));
	assert(b3!=merge(b1, b4));
	assert(b1.merge(b2)==b3);
	assert(b3==merge(merge(b1, b2), b3));
}

void test_binary_hash_bits() {
	// Regression: binary hash longer than 32bits
	binary_hash b=binary_hash("11100110011110000011101110100110001");
	assert(b.test(1));
	assert(b.test(2));
	assert(b.test(3));
	assert(!b.test(4));
	assert(!b.test(5));
}

void test_binary_hash_precision() {
	geolocation l{31.23, 121.473};
	// 1km
	assert(binary_hash_precision(l, 1)==27);
	// 100m
	assert(binary_hash_precision(l, 0.1)==33);
	// 10m
	assert(binary_hash_precision(l, 0.01)==39);
	// 1m
	assert(binary_hash_precision(l, 0.001)==47);
}

void test_binary_encode() {
	geolocation l{31.23, 121.473};
	assert(binary_encode(l, 1)==binary_hash("1"));
	assert(binary_encode(l, 2)==binary_hash("11"));
	assert(binary_encode(l, 3)==binary_hash("111"));
	assert(binary_encode(l, 4)==binary_hash("1110"));
	assert(binary_encode(l, 5)==binary_hash("11100"));
	assert(binary_encode(l, 6)==binary_hash("111001"));
	assert(binary_encode(l, 7)==binary_hash("1110011"));
	assert(binary_encode(l, 8)==binary_hash("11100110"));
	assert(binary_encode(l, 9)==binary_hash("111001100"));
	assert(binary_encode(l, 10)==binary_hash("1110011001"));
	assert(binary_encode(l, 11)==binary_hash("11100110011"));
	assert(binary_encode(l, 12)==binary_hash("111001100111"));
	assert(binary_encode(l, 13)==binary_hash("1110011001111"));
	assert(binary_encode(l, 14)==binary_hash("11100110011110"));
	assert(binary_encode(l, 15)==binary_hash("111001100111100"));
	assert(binary_encode(l, 16)==binary_hash("1110011001111000"));
	assert(binary_encode(l, 17)==binary_hash("11100110011110000"));
	assert(binary_encode(l, 18)==binary_hash("111001100111100000"));
	assert(binary_encode(l, 19)==binary_hash("1110011001111000001"));
	assert(binary_encode(l, 20)==binary_hash("11100110011110000011"));
	assert(binary_encode(l, 21)==binary_hash("111001100111100000111"));
	assert(binary_encode(l, 22)==binary_hash("1110011001111000001111"));
	assert(binary_encode(l, 23)==binary_hash("11100110011110000011110"));
	assert(binary_encode(l, 24)==binary_hash("111001100111100000111100"));
	assert(binary_encode(l, 25)==binary_hash("1110011001111000001111000"));
	assert(binary_encode(l, 26)==binary_hash("11100110011110000011110001"));
	assert(binary_encode(l, 27)==binary_hash("111001100111100000111100010"));
	assert(binary_encode(l, 28)==binary_hash("1110011001111000001111000100"));
	assert(binary_encode(l, 29)==binary_hash("11100110011110000011110001000"));
	assert(binary_encode(l, 30)==binary_hash("111001100111100000111100010001"));
	assert(binary_encode(l, 31)==binary_hash("1110011001111000001111000100011"));
	assert(binary_encode(l, 32)==binary_hash("11100110011110000011110001000110"));
	assert(binary_encode(l, 33)==binary_hash("111001100111100000111100010001100"));
	assert(binary_encode(l, 34)==binary_hash("1110011001111000001111000100011000"));
	assert(binary_encode(l, 35)==binary_hash("11100110011110000011110001000110001"));
	assert(binary_encode(l, 36)==binary_hash("111001100111100000111100010001100011"));
	assert(binary_encode(l, 37)==binary_hash("1110011001111000001111000100011000111"));
	assert(binary_encode(l, 38)==binary_hash("11100110011110000011110001000110001111"));
	assert(binary_encode(l, 39)==binary_hash("111001100111100000111100010001100011111"));
	assert(binary_encode(l, 40)==binary_hash("1110011001111000001111000100011000111111"));
	assert(binary_encode(l, 41)==binary_hash("11100110011110000011110001000110001111111"));
	assert(binary_encode(l, 42)==binary_hash("111001100111100000111100010001100011111111"));
	assert(binary_encode(l, 43)==binary_hash("1110011001111000001111000100011000111111111"));
	assert(binary_encode(l, 44)==binary_hash("11100110011110000011110001000110001111111111"));
	assert(binary_encode(l, 45)==binary_hash("111001100111100000111100010001100011111111110"));
	assert(binary_encode(l, 46)==binary_hash("1110011001111000001111000100011000111111111101"));
	assert(binary_encode(l, 47)==binary_hash("11100110011110000011110001000110001111111111010"));
	assert(binary_encode(l, 48)==binary_hash("111001100111100000111100010001100011111111110100"));
	assert(binary_encode(l, 49)==binary_hash("1110011001111000001111000100011000111111111101000"));
	assert(binary_encode(l, 50)==binary_hash("11100110011110000011110001000110001111111111010001"));
	assert(binary_encode(l, 51)==binary_hash("111001100111100000111100010001100011111111110100010"));
	assert(binary_encode(l, 52)==binary_hash("1110011001111000001111000100011000111111111101000101"));
	assert(binary_encode(l, 53)==binary_hash("11100110011110000011110001000110001111111111010001010"));
	assert(binary_encode(l, 54)==binary_hash("111001100111100000111100010001100011111111110100010101"));
	assert(binary_encode(l, 55)==binary_hash("1110011001111000001111000100011000111111111101000101011"));
	assert(binary_encode(l, 56)==binary_hash("11100110011110000011110001000110001111111111010001010111"));
	assert(binary_encode(l, 57)==binary_hash("111001100111100000111100010001100011111111110100010101111"));
	assert(binary_encode(l, 58)==binary_hash("1110011001111000001111000100011000111111111101000101011111"));
	assert(binary_encode(l, 59)==binary_hash("11100110011110000011110001000110001111111111010001010111110"));
	assert(binary_encode(l, 60)==binary_hash("111001100111100000111100010001100011111111110100010101111100"));
	assert(binary_encode(l, 61)==binary_hash("1110011001111000001111000100011000111111111101000101011111001"));
	assert(binary_encode(l, 62)==binary_hash("11100110011110000011110001000110001111111111010001010111110010"));
	assert(binary_encode(l, 63)==binary_hash("111001100111100000111100010001100011111111110100010101111100101"));
	assert(binary_encode(l, 64)==binary_hash("1110011001111000001111000100011000111111111101000101011111001011"));
}

void test_binary_decode() {
	assert(decode(binary_hash::from_geohash("w")).contains(geolocation{21, 113}));
	assert(decode(binary_hash("11100")).contains(geolocation{21, 113}));
	assert(decode(binary_hash::from_geohash("wt")).contains(geolocation{30.9, 118}));
	assert(decode(binary_hash("1110011001")).contains(geolocation{30.9, 118}));
	assert(decode(binary_hash::from_geohash("wtw")).contains(geolocation{31.6, 121.6}));
	assert(decode(binary_hash("111001100111100")).contains(geolocation{31.6, 121.6}));
	assert(decode(binary_hash::from_geohash("wtw3")).contains(geolocation{31.2, 121.46}));
	assert(decode(binary_hash("11100110011110000011")).contains(geolocation{31.2, 121.46}));
	assert(decode(binary_hash::from_geohash("wtw3r")).contains(geolocation{31.179, 121.619}));
	assert(decode(binary_hash("1110011001111000001110111")).contains(geolocation{31.179, 121.619}));
	assert(decode(binary_hash::from_geohash("wtw3r9")).contains(geolocation{31.1655, 121.624}));
	assert(decode(binary_hash("111001100111100000111011101001")).contains(geolocation{31.1655, 121.624}));
	assert(decode(binary_hash::from_geohash("wtw3r9j")).contains(geolocation{31.1634, 121.6262}));
	assert(decode(binary_hash("11100110011110000011101110100110001")).contains(geolocation{31.1634, 121.6262}));
	assert(decode(binary_hash::from_geohash("wtw3r9jj")).contains(geolocation{31.16366, 121.62569}));
	assert(decode(binary_hash("1110011001111000001110111010011000110001")).contains(geolocation{31.16366, 121.62569}));
	assert(decode(binary_hash::from_geohash("wtw3r9jjz")).contains(geolocation{31.163728, 121.625841}));
	assert(decode(binary_hash("111001100111100000111011101001100011000111111")).contains(geolocation{31.163728, 121.625841}));
	assert(decode(binary_hash::from_geohash("wtw3r9jjzy")).contains(geolocation{31.1637416, 121.625857}));
	assert(decode(binary_hash("11100110011110000011101110100110001100011111111110")).contains(geolocation{31.1637416, 121.625857}));
	assert(decode(binary_hash::from_geohash("wtw3r9jjzyj")).contains(geolocation{31.1637396, 121.6258588}));
	assert(decode(binary_hash("1110011001111000001110111010011000110001111111111010001")).contains(geolocation{31.1637396, 121.6258588}));
	assert(decode(binary_hash::from_geohash("wtw3r9jjzyjc")).contains(geolocation{31.16373922, 121.62585927}));
	assert(decode(binary_hash("111001100111100000111011101001100011000111111111101000101011")).contains(geolocation{31.16373922, 121.62585927}));
	assert(!decode(binary_hash::from_geohash("wtw3r9jjzyjc")).contains(geolocation{31.16374922, 121.62585927}));
}

void test_binary_neighbor() {
	binary_hash b("11100110");
	assert(neighbor(b, {-1, -1})==binary_hash("11100001"));
	assert(neighbor(b, {-1,  0})==binary_hash("11100011"));
	assert(neighbor(b, {-1,  1})==binary_hash("11101001"));
	assert(neighbor(b, { 0, -1})==binary_hash("11100100"));
	assert(neighbor(b, { 0,  1})==binary_hash("11101100"));
	assert(neighbor(b, { 1, -1})==binary_hash("11100101"));
	assert(neighbor(b, { 1,  0})==binary_hash("11100111"));
	assert(neighbor(b, { 1,  1})==binary_hash("11101101"));
}

void test_encode() {
	geolocation l{31.16373922, 121.62585927};
	assert(encode(l, 1)=="w");
	assert(encode(l, 2)=="wt");
	assert(encode(l, 3)=="wtw");
	assert(encode(l, 4)=="wtw3");
	assert(encode(l, 5)=="wtw3r");
	assert(encode(l, 6)=="wtw3r9");
	assert(encode(l, 7)=="wtw3r9j");
	assert(encode(l, 8)=="wtw3r9jj");
	assert(encode(l, 9)=="wtw3r9jjz");
	assert(encode(l, 10)=="wtw3r9jjzy");
	assert(encode(l, 11)=="wtw3r9jjzyj");
	assert(encode(l, 12)=="wtw3r9jjzyjc");
}

void test_decode() {
	assert(decode("w").contains(geolocation{21, 113}));
	assert(decode("wt").contains(geolocation{30.9, 118}));
	assert(decode("wtw").contains(geolocation{31.6, 121.6}));
	assert(decode("wtw3").contains(geolocation{31.2, 121.46}));
	assert(decode("wtw3r").contains(geolocation{31.179, 121.619}));
	assert(decode("wtw3r9").contains(geolocation{31.1655, 121.624}));
	assert(decode("wtw3r9j").contains(geolocation{31.1634, 121.6262}));
	assert(decode("wtw3r9jj").contains(geolocation{31.16366, 121.62569}));
	assert(decode("wtw3r9jjz").contains(geolocation{31.163728, 121.625841}));
	assert(decode("wtw3r9jjzy").contains(geolocation{31.1637416, 121.625857}));
	assert(decode("wtw3r9jjzyj").contains(geolocation{31.1637396, 121.6258588}));
	assert(decode("wtw3r9jjzyjc").contains(geolocation{31.16373922, 121.62585927}));
	assert(!decode("wtw3r9jjzyjc").contains(geolocation{31.16374922, 121.62585927}));
}

void test_encode_precision_range() {
	typedef std::vector<std::string> hs;
	geolocation l{31.16373922, 121.62585927};
	hs codes;
	encode_precision_range(l, codes, 1, 9);
	assert(codes==hs({"w", "wt", "wtw", "wtw3", "wtw3r", "wtw3r9", "wtw3r9j", "wtw3r9jj", "wtw3r9jjz"}));
}

void test_hash_precision() {
	geolocation l{31.23, 121.473};
	// 1km
	assert(hash_precision(l, 1)==5);
	// 100m
	assert(hash_precision(l, 0.1)==6);
	// 10m
	assert(hash_precision(l, 0.01)==7);
	// 1m
	assert(hash_precision(l, 0.001)==9);
}

void test_base_hash() {
	geolocation l{31.23, 121.473};
	// 1km
	assert(base_hash(l, 1)=="wtw3s");
	// 100m
	assert(base_hash(l, 0.1)=="wtw3sj");
	// 10m
	assert(base_hash(l, 0.01)=="wtw3sjj");
	// 1m
	assert(base_hash(l, 0.001)=="wtw3sjjzy");
}

void test_hash_contains() {
	assert(hash_contains("w", geolocation{21, 113}));
	assert(hash_contains("wt", geolocation{30.9, 118}));
	assert(hash_contains("wtw", geolocation{31.6, 121.6}));
	assert(hash_contains("wtw3", geolocation{31.2, 121.46}));
	assert(hash_contains("wtw3r", geolocation{31.179, 121.619}));
	assert(hash_contains("wtw3r9", geolocation{31.1655, 121.624}));
	assert(hash_contains("wtw3r9j", geolocation{31.1634, 121.6262}));
	assert(hash_contains("wtw3r9jj", geolocation{31.16366, 121.62569}));
	assert(hash_contains("wtw3r9jjz", geolocation{31.163728, 121.625841}));
	assert(hash_contains("wtw3r9jjzy", geolocation{31.1637416, 121.625857}));
	assert(hash_contains("wtw3r9jjzyj", geolocation{31.1637396, 121.6258588}));
	assert(hash_contains("wtw3r9jjzyjc", geolocation{31.16373922, 121.62585927}));
	assert(!hash_contains("wtw3r9jjzyjc", geolocation{31.16374922, 121.62585927}));
}

void test_neighbor() {
	assert(neighbor("wtw3s", {-1, -1})=="wtw37");
	assert(neighbor("wtw3s", {-1, 0})=="wtw3k");
	assert(neighbor("wtw3s", {-1, 1})=="wtw3m");
	assert(neighbor("wtw3s", {0, -1})=="wtw3e");
	assert(neighbor("wtw3s", {0, 1})=="wtw3t");
	assert(neighbor("wtw3s", {1, -1})=="wtw3g");
	assert(neighbor("wtw3s", {1, 0})=="wtw3u");
	assert(neighbor("wtw3s", {1, 1})=="wtw3v");
	assert(neighbor("wtw3sjj", {-1, -1})=="wtw3shu");
	assert(neighbor("wtw3sjj", {-1, 0})=="wtw3shv");
	assert(neighbor("wtw3sjj", {-1, 1})=="wtw3shy");
	assert(neighbor("wtw3sjj", {0, -1})=="wtw3sjh");
	assert(neighbor("wtw3sjj", {0, 1})=="wtw3sjn");
	assert(neighbor("wtw3sjj", {1, -1})=="wtw3sjk");
	assert(neighbor("wtw3sjj", {1, 0})=="wtw3sjm");
	assert(neighbor("wtw3sjj", {1, 1})=="wtw3sjq");
}

void test_hash_codes() {
	typedef std::vector<std::string> hs;
	geolocation l{31.23, 121.473};
	{
		hs codes;
		hash_codes(l, 1, codes);
		assert(codes==hs({"wtw37", "wtw3k", "wtw3m", "wtw3e", "wtw3t", "wtw3g", "wtw3u", "wtw3v", "wtw3s"}));
	}
	{
		hs codes;
		hash_codes(l, 0.1, codes);
		assert(codes==hs({"wtw3eu", "wtw3sh", "wtw3sk", "wtw3ev", "wtw3sm", "wtw3ey", "wtw3sn", "wtw3sq", "wtw3sj"}));
	}
	{
		hs codes;
		hash_codes(l, 0.01, codes);
		assert(codes==hs({"wtw3shu", "wtw3shv", "wtw3shy", "wtw3sjh", "wtw3sjn", "wtw3sjk", "wtw3sjm", "wtw3sjq", "wtw3sjj"}));
	}
	{
		hs codes;
		hash_codes(l, 0.001, codes);
		assert(codes==hs({"wtw3sjjzt", "wtw3sjjzw", "wtw3sjjzx", "wtw3sjjzv", "wtw3sjjzz", "wtw3sjmbj", "wtw3sjmbn", "wtw3sjmbp", "wtw3sjjzy"}));
	}
}

int main() {
	test_geolocation();
	test_bbox1();
	test_bbox2();
	test_binary_hash_bits();
	test_binary_hash_precision();
	test_binary_encode();
	test_binary_decode();
	test_binary_neighbor();
	test_encode();
	test_decode();
	test_encode_precision_range();
	test_hash_precision();
	test_base_hash();
	test_hash_contains();
	test_neighbor();
	test_hash_codes();
	return 0;
}