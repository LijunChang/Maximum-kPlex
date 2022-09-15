#ifndef _ARRAY_H
#define _ARRAY_H
#endif // _ARRAY_H

#include <cstring>

class Array
{
protected:
    int *array;
    int *index_in_array;
    int array_size;
    int array_capacity;

public:
    Array()
    {
        array = NULL;
        index_in_array = NULL;
        array_size = 0;
        array_capacity = 0;
    }

    Array(int capacity)
    {
        array_capacity = capacity + 1;
        array = new int[array_capacity];
        index_in_array = new int[array_capacity];
        //memset(array, 0, sizeof(int) * (array_capacity));
        //memset(index_in_array, 0, sizeof(int) * (array_capacity));
        for(int i = 0; i < array_capacity; i++)
        {
            array[i] = -1;
            index_in_array[i] = -1;
        }
        array_size = 0;
    }

    ~Array()
    {
        delete[] array;
        array = NULL;
        delete[] index_in_array;
        index_in_array = NULL;
        array_size = 0;
        array_capacity = 0;
    }

    void clear()
    {
        for(int i = 0; i < array_size; i++)
            index_in_array[array[i]] = -1;
        //for(int i = 0; i < array_capacity; i++)
            //index_in_array[i] = 0;
        array_size = 0;
    }

    void insert_element(int e)
    {
        array[array_size] = e;
        index_in_array[e] = array_size;
        array_size++;
    }

    void delete_element(int e)
    {
        array_size--;
        int last_e = array[array_size];
        int e_pos = index_in_array[e];
        array[e_pos] = last_e;
        index_in_array[last_e] = e_pos;
        index_in_array[e] = -1;
    }

    int index_of(int e)
    {
        return index_in_array[e];
    }

    int element_at(int i)
    {
        return array[i];
    }

    bool is_in_array(int e)
    {
        if(index_in_array[e] == -1)
            return false;
        else
            return true;
    }

    int size()
    {
        return array_size;
    }

    int begin()
    {
        return 0;
    }

    int end()
    {
        return array_size;
    }

    int rand_element()
    {
        return array[rand() % array_size];
    }
};
