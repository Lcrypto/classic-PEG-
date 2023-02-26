#include "ACE_toolbox.h"

bool ACE_detection(std::vector<std::vector<int>> vn_neighbor, std::vector<std::vector<int>> cn_neighbor, int d_ace, int eta_ace, int starting_v)
{
    /* 
    This algorithm analyze the ACE peropety given a matrix
    One has to provide the cn indecs (starting from 0) conncoted to each vn
    One has to provide the vn indeces (staring from 0) connected to each cn
    One has to prvoide two critical facors d_ace and eta_ace
    One has to provide the variable node it exmines, see paper in detail if you want to know the exat meaning of these two dectors 
    Function returns true, if the metrix passes ace dectction,
    Function returns false, otherwise
    */

    //Initialization
    std::vector<std::vector<int>> node_ace(2); //[0] -> check node ; [1] -> variable node
    std::vector<std::vector<int>> p(2);        //[0] -> check node; [1] -> variable node
    int p_tempt;
    int zero = 0;
    std::vector<std::vector<int>> node_set_kids, node_set_parents;
    std::vector<int> his_kids, this_kid;
    bool not_appear;
    int node_type_kid, node_type_parent;

    //Do some precomputation
    node_ace[0].assign(cn_neighbor.size(), 0);
    node_ace[1].assign(vn_neighbor.size(), 0);
    for (unsigned ii = 0; ii < vn_neighbor.size(); ii++)
    {
        node_ace[1][ii] = std::max(zero, (int)vn_neighbor[ii].size() - 2);
    }


    // display(node_ace);
    // std::cout<<"-------"<<std::endl;
    p[0].assign(cn_neighbor.size(), 10000);
    p[1].assign(vn_neighbor.size(), 10000);
    p[1][starting_v]=node_ace[1][starting_v];
    node_set_parents.clear();
    node_set_parents.push_back({-1,starting_v});
    for (int l = 1; l <= d_ace; l++)
    {
        if (l % 2 == 0)
        {
            node_type_parent = 0;
            node_type_kid = 1;
        }
        else
        {
            node_type_parent = 1;
            node_type_kid = 0;
        }
        if (node_set_parents.size() == 0)
        {
            std::cout << "There is no node in parent set. Returned Success" << std::endl;
            return true;
        }
        
        else
        {
            node_set_kids.clear();
            for (unsigned kk=0;kk<node_set_parents.size();kk++)
            {
                if (node_type_parent == 0)
                {
                    his_kids = cn_neighbor[node_set_parents[kk][1]];
                }
                else
                {
                    his_kids = vn_neighbor[node_set_parents[kk][1]];
                }
                for (const auto this_kid : his_kids)
                {
                    if (this_kid != node_set_parents[kk][0])
                    {
                        
                        p_tempt = p[node_type_parent][node_set_parents[kk][1]] + node_ace[node_type_kid][this_kid];
                        if (p_tempt + p[node_type_kid][this_kid] - node_ace[1][starting_v] - node_ace[node_type_kid][this_kid] < eta_ace)
                        {
                            std::cout << "Info: find ACE=" << p_tempt + p[node_type_kid][this_kid] - node_ace[1][starting_v] - node_ace[node_type_kid][this_kid] << ". Exit.." << std::endl;
                            return false;
                        }
                        else
                        {
                            if (p_tempt < p[node_type_kid][this_kid])
                            {
                                p[node_type_kid][this_kid]=p_tempt;
                                not_appear = true;
                                for (unsigned jj = 0; jj < node_set_kids.size(); jj++)
                                {
                                    if (node_set_kids[jj][1] == this_kid)
                                    {
                                        not_appear = false;
                                        break;
                                    }
                                }
                                if (not_appear)
                                {
                                    node_set_kids.push_back({node_set_parents[kk][1],this_kid});
                                }
                            }
                        }
                    }
                }
            }
        }
        //display part
        // std::cout<<"l="<<l<<". "<<std::endl;
        // std::cout<<"kid set is :"<< std::endl;
        // for(const auto aa: node_set_kids)
        //     std::cout<<aa[1]<<" ";
        // std::cout<<std::endl;
        // std::cout<<"p value is :"<<std::endl;
        // for(const auto aa: node_set_kids)
        //     std::cout<<p[node_type_kid][aa[1]]<<" ";
        // std::cout<<std::endl;
        node_set_parents=node_set_kids;
    }
    //std::cout<<"Info: ACE dectection passed ..."<<starting_v<<std::endl;
    return true;
}


std::vector<std::vector<int>> find_cn_neighbors(std::vector<std::vector<int>> & input_matrix)
{
    /*
    This function finds all variable nodes connect to each check node
    if one check node currently does not conncet variable node
    we will also keep a position to this check node 
    but vector size will be empet, i.e., 0 
    */
   std::vector<std::vector<int>> cn_neighbor(input_matrix.size());
   unsigned vari_num = input_matrix[0].size();
   for (unsigned ii=0; ii< cn_neighbor.size() ; ii++)
   {
       cn_neighbor[ii].clear();
       for (unsigned jj=0;jj<vari_num;jj++)
       {
           if(input_matrix[ii][jj]!=0)
           {
               cn_neighbor[ii].push_back(jj);
           }
       }
   }
   return cn_neighbor;
}

std::vector<std::vector<int>> find_vn_neighbers(std::vector<std::vector<int>> & input_matrix)
{
    /*
    find all check nodee connets to each variable node
    rule is same to find_cn_neighbors 
    */
   unsigned vari_num = input_matrix[0].size();
   unsigned check_num = input_matrix.size();
   std::vector<std::vector<int>> vn_neighbers(vari_num);
   for (unsigned ii=0 ;ii< vari_num;ii++)
   {
       vn_neighbers[ii].clear();
       for (unsigned jj=0;jj<check_num;jj++)
       {
           if(input_matrix[jj][ii]!=0)
           {
               vn_neighbers[ii].push_back(jj);
           }
       }
   }
   return vn_neighbers;

}

int rankOfMatrix(std::vector<std::vector<int>> mat) 
{ 
    int R = (int) mat.size();
    int C = (int) mat[0].size();

    int rank = C; 
    //std::cout<<R<<" "<<C<<std::endl;
    for (int row = 0; row < rank; row++) 
    { 
        
        if(row==R)
        {
            
            return R;
        }
        
        // Before we visit current row 'row', we make 
        // sure that mat[row][0],....mat[row][row-1] 
        // are 0. 

        // Diagonal element is not zero
        if (mat[row][row])
        {
            for (int col = 0; col < R; col++)
            {
                if (col != row)
                {
                    // This makes all entries of current
                    // column as 0 except entry 'mat[row][row]'
                    // column as 0 except entry 'mat[row][row]'
                    if (mat[col][row]==1)
                    {
                        for( int i=row ; i<rank ;i ++ )
                        {
                            mat[col][i] = (mat[col][i]+mat[row][i])%2;
                        }
                    }
                }
            }
        }

        // Diagonal element is already zero. Two cases 
        // arise: 
        // 1) If there is a row below it with non-zero 
        //    entry, then swap this row with that row 
        //    and process that row 
        // 2) If all elements in current column below 
        //    mat[r][row] are 0, then remvoe this column 
        //    by swapping it with last column and 
        //    reducing number of columns by 1. 
        else
        { 
            bool reduce = true; 
  
            /* Find the non-zero element in current 
                column  */
            for (int i = row + 1; i < R;  i++) 
            { 
                // Swap the row with non-zero element 
                // with this row. 
                if (mat[i][row]) 
                { 
                    swap(mat, row, i, rank); 
                    reduce = false; 
                    //std::cout<<"here"<<std::endl;
                    break ; 
                } 
            }  
            // If we did not find any row with non-zero 
            // element in current columnm, then all 
            // values in this column are 0. 
            if (reduce) 
            { 
                // Reduce number of columns 
                rank--;  
                // Copy the last column here 
                for (int i = 0; i < R; i ++) 
                    mat[i][row] = mat[i][rank]; 
            }
            // Process this row again
            row--;
        }

        // Uncomment these lines to see intermediate results
        //  display(mat, R, C);
        //  printf("\n");
        //std::cout<<row<<std::endl;
        // if (C==1032)
        // {
        //     std::cout<<row<<std::endl;
        // }
        // if(row==257&&C==1032)
        // {
        //     std::ofstream myfile("check.txt");
        //     for(auto aa:mat)
        //     {
        //         for(auto bb:aa)
        //         {
        //             myfile<<bb<<"  ";
        //         }
        //         myfile<<std::endl;
        //     }
        //     myfile.close();
        // }
    }
    return rank;
}

void swap(std::vector<std::vector<int>> & mat, int row1, int row2,  int col) 
{ 
    for (int i = 0; i < col; i++) 
    { 
        int temp = mat[row1][i]; 
        mat[row1][i] = mat[row2][i]; 
        mat[row2][i] = temp; 
    } 
} 

void display(std::vector<std::vector<int>> mat, int row, int col) 
{ 
    for (int i = 0; i < row; i++) 
    { 
        for (int j = 0; j < col; j++) 
            printf("  %d", mat[i][j]); 
        printf("\n"); 
    } 
} 

void display(std::vector<std::vector<int>> mat)
{
    for (unsigned i = 0; i < mat.size(); i++)
    {
        if (mat[i].size() == 0)
        {
            std::cout << "*****empty line*****" << std::endl;
        }
        else
        {
            for (unsigned j = 0; j < mat[i].size(); j++)
                printf("  %d", mat[i][j]);
            printf("\n");
        }
    }
}

void display(std::vector<int> mat)
{
    for(const auto aa:mat)
        std::cout<<aa <<"  ";
    std::cout<<std::endl;
}