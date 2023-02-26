# include<read_H_matrix.h>


void get_whole_H(std::vector<std::vector<int>>&  H,std::string filename, const char type[])
{
    // this function get M*N whold parity check matrix
    // it can get H from different tyope of files
    if(strcmp(type,"classic")==0)
    {
        std::ifstream myfile(filename);
        if(myfile.is_open())
        {
            int row, column;
            myfile>>column>>row;
            H.assign(row, std::vector<int> (column,0));
            for (int ii = 0; ii<row ;ii++)
            {
                for(int jj=0;jj<column;jj++)
                {
                    myfile>>H[ii][jj];
                }
            }
            std::cout<<"File :" << filename<<" is successfully read." <<std::endl;           
        }
        else
        {
            std::cout<<"can't open file :" << filename<<". check again plz." <<std::endl;
        }
        
    }
    else
    {
        std::cout<<"Info: No such type: "<<type<<". Check again, please. "<<std::endl;
    }   
}


void write_matrix(std::vector<std::vector<int>> &H,std::string filename, const char type[], int cir_size)
{
    // this function get M*N whold parity check matrix
    // it can write H to a text file with different type of files
    if (strcmp(type, "circulant")==0)
    {        
        /*  line 1 :    circulant (lifting) size, the number of horizontal layers, the number of  vertical layers
            line 2-end: row-by-row, the first row of each circulant*/
        int row_layer_num = (int)H.size()/cir_size;
        int col_layer_num = (int)H[0].size()/cir_size;
        std::ofstream filehandle(filename);
        filehandle<<cir_size<<" "<<row_layer_num<<" "<<col_layer_num<<std::endl;
        for(int row_ind=0; row_ind<row_layer_num; row_ind++)
        {
            for(int col_ind=0; col_ind<col_layer_num; col_ind++)
            {
                int start_row_ind = row_ind*cir_size;
                int start_col_ind = col_ind*cir_size;

                for (int ii=0; ii<cir_size; ii++)
                {
                    filehandle<<H[start_row_ind][start_col_ind+ii]<<" ";
                }

                filehandle<<std::endl;
            }
        }
        filehandle.close();
    }
    else if (strcmp(type, "protomatrix")==0)
    {
        // Assume no parallel edges!  
        int row_layer_num = (int)H.size()/cir_size;
        int col_layer_num = (int)H[0].size()/cir_size;
        std::ofstream filehandle(filename);
        filehandle<<cir_size<<" "<<row_layer_num<<" "<<col_layer_num<<std::endl;
        for(int row_ind=0; row_ind<row_layer_num; row_ind++)
        {
            for(int col_ind=0; col_ind<col_layer_num; col_ind++)
            {
                int start_row_ind = row_ind*cir_size;
                int start_col_ind = col_ind*cir_size;

                for (int ii=0; ii<cir_size; ii++)
                {
                    if(H[start_row_ind][start_col_ind+ii]==1)
                    {
                        filehandle<<ii<<" ";
                        break;
                    }
                    if (ii==cir_size-1)
                    {
                        filehandle<<-1<<" ";
                    }
                }
                
            }
            filehandle<<std::endl;
        }
        filehandle.close();
    }
    else if (strcmp(type,"chinn")==0||strcmp(type,"Chinn")==0)
    {
        std::cout<<"The author is too lazy to write this part ... "<<std::endl;
        std::cout<<"Chinn has the follwing format ..." <<std::endl;
        std::cout<<"Line 1: Number of variable nodes.\n";
        std::cout<<"Line 2: Line 2: Number of check nodes.\n";
        std::cout<<"Line 3: Number of edges in the graph.\n";
        std::cout<<"Line 4: Field size.\n";
        std::cout<<"Line 5: Max degree of variable nodes.\n";
        std::cout<<"Line 6: Max degree of check nodes.\n";
        std::cout<<"Line 7 through end: <vnode index><space><cnode index>\n";
    }
    else
    {
        std::cout<<"From function write_matrix: No such type: "<<type<<". Check again, please. "<<std::endl;
    }

}