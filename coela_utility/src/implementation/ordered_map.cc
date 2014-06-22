#include "../ordered_map.h"


//#include "../level1/my_debug_config.h"

using namespace std;
#include <stdexcept>

namespace coela {



//=================================================================================
//Ordered map:

void OrderedMap::clear()
{
    table.clear();
    index_map.clear();
}

void OrderedMap::add_key(const KeyValueComment& input)
{
    if (key_exists(input.key) && !input.key.empty()) {
      throw runtime_error(
          "Ordered map: Cannot add a key that already exists (key:"+input.key+")");
    }
    else {
        table.push_back(input);
        index_map[input.key]=table.size()-1;
    }
}

void OrderedMap::remove_key(const string& key)
{
    if (key_exists(key) && !key.empty()) {
        //comments have empty key, would be impossible to know which one to erase.
        table.erase(table.begin()+index_map[key]);
        index_map.erase(key);
        update_index_map();
    }
}


bool OrderedMap::key_exists(const string& key) const
{
    if (index_map.find(key)!=index_map.end()) { return true; }
    else { return false; }
}

void OrderedMap::merge_with_prefix_map(const OrderedMap& prefix_map)
{
    OrderedMap m(prefix_map);
    for (vector<KeyValueComment>::const_iterator i=table.begin();
            i!=table.end(); ++i) {
        if (!m.key_exists(i->key)) {
            m.add_key(*i);
        }
    }
    *this=m;
}


const string& OrderedMap::operator[](const string& key) const
{
    if (key_exists(key)) { return table[ index_map.find(key)->second ].value; }
    throw domain_error("Ordered map op[]- Key: "+key +" does not exist");
}
string & OrderedMap::operator[](const string& key)
{
    if (key_exists(key) ==false) { add_key(KeyValueComment(key,"")); }
    return table[index_map[key]].value;
}

void OrderedMap::update_index_map()
{
    for (size_t i=0; i!=table.size(); ++i) {
        index_map[table[i].key]=i;
    }
}

//End ordered map
//=================================================================================

}//end namespace coela
