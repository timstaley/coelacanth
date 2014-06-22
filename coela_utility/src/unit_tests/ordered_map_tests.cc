/*
 * File:   FitsHeader_tests.cc
 * Author: ts337
 *
 * Created on 07 February 2011, 12:41
 */
#include <UnitTest++/UnitTest++.h>
#include <stdexcept>
#include "../ordered_map.h"
#include <iostream>

using namespace coela;
using namespace std;

SUITE(Ordered_Map)
{
    TEST(Notify_Suite_Has_Been_Run) {
        cout << "*** \"OrderedMap\" unit tests running ***" <<endl;
    }

    TEST(Basic_Functionality) {
        OrderedMap map1;
        key_value_comment k0("somekey","someval");
        key_value_comment k1("anotherkey","someval");
        key_value_comment k2("key3","value3");

        vector<string> value_strings;

        map1.add_key(k0);
        value_strings.push_back(k0.value);

        map1.add_key(k1);
        value_strings.push_back(k1.value);

        map1.add_key(k2);
        value_strings.push_back(k2.value);


        CHECK_THROW(map1.add_key(k0), std::runtime_error); //Cannot add 2 entries with same key

        CHECK_EQUAL(value_strings.size(),map1.size());

        for (size_t index=0; index!=map1.size(); ++index) {
            CHECK_EQUAL(value_strings[index], map1[index].value);
        }

        map1.remove_key(k1.key);
        CHECK_EQUAL(value_strings.size()-1,map1.size());

        OrderedMap map2;
        map2.add_key(k0);
        map2.add_key(k2);
        CHECK(map1 == map2);
    }

    TEST(Merge_with_prefix_map) {
        OrderedMap map1;
        key_value_comment k0("somekey","someval");
        key_value_comment k1("anotherkey","someval");
        key_value_comment k2("key3","value3");

        map1.add_key(k0);
        map1.add_key(k1);
        map1.add_key(k2);

        key_value_comment k3("key4","value4");

        OrderedMap prefix_map;
        prefix_map.add_key(k1);
        prefix_map.add_key(k3);

        OrderedMap manual_merge;
        manual_merge.add_key(k1);
        manual_merge.add_key(k3);

        manual_merge.add_key(k0);
        manual_merge.add_key(k2);

        map1.merge_with_prefix_map(prefix_map);

        CHECK_EQUAL(manual_merge.size(), map1.size());
        CHECK(map1 == manual_merge);


    }
}


