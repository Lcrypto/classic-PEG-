#include "CPEG_toolbox.h"

std::vector<std::vector<int>> single_circulant_generator(int num, int cirsize)
{
    std::vector<std::vector<int>> circulant(cirsize, std::vector<int>(cirsize,0));
    std::random_device device;
    std::mt19937 generator(device());
    std::uniform_int_distribution<int> distribution(0,cirsize-1);

    std::vector<int> random_set;
    int cur_ran;
    bool exist=false;
    bool not_finish=true;
    int cur_point;
    do
    {
        cur_ran=distribution(generator);
        exist=false;
        for (unsigned ii=0; ii<random_set.size();ii++)
        {
            if(random_set[ii]==cur_ran)
            {
                exist=true;
                break;
            }
        }
        if(!exist)
        {
            random_set.push_back(cur_ran);
        }
        if((int)random_set.size()==num)
        {
            not_finish=false;
        }

    } while(not_finish);

    for (int ii=0; ii<num; ii++)
    {
        cur_point=random_set[ii];
        for (int jj=0; jj<cirsize; jj++)
        {
            circulant[jj][cur_point]=1;
            cur_point=(cur_point+1)%cirsize;
        }
    }

    return circulant;
}

std::vector<std::vector<int>> identity_circulant_generator(int cirsize)
{
    std::vector<std::vector<int>> circulant(cirsize,std::vector<int>(cirsize,0));
    for (int ii=0; ii<cirsize;ii++)
    {
        circulant[ii][ii]=1;
    }
    return circulant;
}

std::vector<std::vector<int>> column_circulant_generator(std::vector<int> column_vector, int cirsize)
{
    int start,end;
    int row;
    std::vector<std::vector<int>> cur_circulant;
    std::vector<std::vector<int>> column_circulant((int)column_vector.size()*cirsize, std::vector<int>(cirsize,0));
    for (unsigned ii=0 ; ii<column_vector.size();ii++)
    {
        start=ii*cirsize;
        end=(ii+1)*cirsize-1;
        if(column_vector[ii]>0)
        {
            cur_circulant=single_circulant_generator(column_vector[ii], cirsize);
            row=0;
            for( int jj=start; jj<=end; jj++)
            {
                column_circulant[jj]=cur_circulant[row];
                row++;

            }
        }
    }
    return column_circulant;
}

std::vector<std::vector<int>> column_circulant_generator(std::vector<int> column_vector, int cirsize, const char type[])
{
    int start, end;
    int row;
    std::vector<std::vector<int>> cur_circulant;
    std::vector<std::vector<int>> column_circulant((int)column_vector.size() * cirsize, std::vector<int>(cirsize, 0));
    if (strcmp(type, "identity") == 0)
    {
        for (unsigned ii = 0; ii < column_vector.size(); ii++)
        {
            start = ii * cirsize;
            end = (ii+1) * cirsize - 1;
            if (column_vector[ii] > 0)
            {
                if(column_vector[ii]!=1)
                {
                    std::cout<<"Info: Illegle use a non_zero element to generate a identiry ..."<<std::endl;
                }
                cur_circulant = identity_circulant_generator(cirsize);
                row = 0;
                for (int jj = start; jj <= end; jj++)
                {
                    column_circulant[jj] = cur_circulant[row];
                    row++;
                }
            }
        }
    }
    else
    {
        std::cout<<"Info: No such type, please check again ..." <<std::endl;
    }
    return column_circulant;
    
}

void add_new_colum_circulant(std::vector<std::vector<int>> & parity_check, std::vector<std::vector<int>> &new_column_circulant)
{
    for(unsigned ii=0 ; ii<parity_check.size(); ii++)
    {
        std::copy(new_column_circulant[ii].begin(),new_column_circulant[ii].end(),std::back_inserter(parity_check[ii]));
    }
}