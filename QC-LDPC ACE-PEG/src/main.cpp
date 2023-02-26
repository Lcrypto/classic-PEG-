using namespace std; 
#include <ACE_toolbox.h>
#include <CPEG_toolbox.h>
#include "PBRL_LDPC_H_Generator.h"
#include "read_H_matrix.h"

// Driver program to test above functions
int main(int argc, char *argv[])
{

    std::string description_filename = "none.txt";
    std::string output_filename = "none.txt";
    std::string output_type = "protomatrix";

    int high_rate_row_ind=0;
    int high_rate_column_ind=3;
    int d_ace=3;
    int eta_ace=2;
    int pre_lifting_size=8;
    int cir_size= 32;
    int num_of_in = (argc - 1) / 2;
    std::vector<std::vector<int>> proto_matrix;
    for (int i = 0; i < num_of_in; ++i)
    {
        if (strcmp(argv[2 * i + 1], "-description_file") == 0)
        {
            description_filename = argv[2 * i + 2];
        }
        else if (strcmp(argv[2 * i + 1], "-output_file") == 0)
        {
            output_filename = argv[2 * i + 2];
        }
        else if (strcmp(argv[2 * i + 1], "-output_type") == 0)
        {
            output_type = argv[2 * i + 2];
        }
        else
        {
            std::cout << argv[2 * i + 1] << " is not recgonized ... plz check again" << std::endl;
            return 0;
        }
    }

    // This part reads the description file
    std::ifstream filehandle(description_filename);
    if(filehandle.is_open())
    {
        // read the matrix info .. 
        int row_num, col_num;
        filehandle>>row_num;
        filehandle>>col_num;
        proto_matrix.assign(row_num, std::vector<int>(col_num,0));
        for (int row_ind = 0; row_ind<row_num; row_ind++)
        {
            for(int col_ind = 0; col_ind<col_num; col_ind++)
            {
                filehandle>>proto_matrix[row_ind][col_ind];
            }
        }
        filehandle >> high_rate_row_ind;
        filehandle >> high_rate_column_ind;
        filehandle >> d_ace;
        filehandle >> eta_ace;
        filehandle >> pre_lifting_size;
        filehandle >> cir_size;
        std::cout << "The description file: " << description_filename << " has been loaded successfully" << std::endl;
        std::vector<std::vector<int>> pre_lifted_matrix = pre_lifter(proto_matrix, pre_lifting_size, high_rate_column_ind);
        std::vector<std::vector<int>> matrix = ACE_PEG_generator(pre_lifted_matrix, (high_rate_row_ind + 1) * pre_lifting_size - 1, (high_rate_column_ind + 1) * pre_lifting_size - 1, cir_size, d_ace, eta_ace);
        write_matrix(matrix,output_filename,output_type.data(),cir_size);
    }
    else
    {
        std::cout<<"The description file: "<<description_filename<<" can not be open."<<std::endl;
    }





    return 0;
}

// int main()
// {

//     // std::vector<std::vector<int>> protomatrix{
//     //     {3, 3, 3, 3, 1, 1, 1, 2, 1, 1, 0, 0, 0},
//     //     {1, 1, 1, 1, 3, 3, 3, 1, 2, 2, 0, 0, 0},
//     //     {2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0},
//     //     {1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0},
//     //     {2, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1}
//     //     };

//     std::vector<std::vector<int>> protomatrix{{7, 3, 2, 2}};

    
//     int high_rate_row_ind=0;
//     int high_rate_column_ind=3;
//     int d_ace=3;
//     int eta_ace=2;
//     int pre_lifting_size=8;
//     int cir_size= 32;

//     std::vector<std::vector<int>> pre_lifted_matrix = pre_lifter(protomatrix,pre_lifting_size,high_rate_column_ind);
//     std::vector<std::vector<int>> maxrix=ACE_PEG_generator(pre_lifted_matrix,(high_rate_row_ind+1)*pre_lifting_size-1,(high_rate_column_ind+1)*pre_lifting_size-1,cir_size,d_ace,eta_ace);

//     // rank of matrix
//     // std::vector<std::vector<int>> H;
//     // int rank;
//     // std::string filename="H_minLUT_classic.txt";
//     // get_whole_H(H,filename,"classic");
//     // rank = rankOfMatrix(H);
//     // std::cout<<"rank: "<<rank<<std::endl;

//     // std::vector<std::vector<int>> mat{
//     //     {1,1,0,0},
//     //     {1,1,0,1},
//     //     {0,0,1,1},
//     //     {1,1,1,1}
//     // };
//     // int a=rankOfMatrix(mat);
//     return 0;
// }

// // std::cout << "test generating random desigend circulant.." << std::endl;
// // std::vector<std::vector<int>> test1 = single_circulant_generator(2, 8);
// // std::cout << "result is:" << std::endl;
// // display(test1);
// // std::cout << "test idenity matrix generator.." << std::endl;
// // std::vector<std::vector<int>> test2 = identity_circulant_generator(8);
// // std::cout << "result is:" << std::endl;
// // display(test2);

// // std::vector<int> column_vector1{0,1,2};
// // std::vector<int> column_vector2{0,1,0};
// // std::cout<<"now testing random column circulant generator ..." <<std::endl;
// // std::vector<std::vector<int>> column_circulant1=column_circulant_generator(column_vector1,6);
// // display(column_circulant1);
// // std::cout<<"--------------------------------"<<std::endl;
// // std::cout<<"now testing identity column circulant generator ..." <<std::endl;
// // std::vector<std::vector<int>> column_circulant2=column_circulant_generator(column_vector2,6,"identity");
// // display(column_circulant2);
// // std::cout<<"--------------------------------"<<std::endl;
// // std::cout<<"now testing identity column circulant combination ..." <<std::endl;
// // add_new_colum_circulant(column_circulant1,column_circulant2);
// // display(column_circulant1);

// // std::cout << "the matrix we used is..." << std::endl;
// // display(matrix);
// // std::cout << "-------------------------" << std::endl;
// // std::cout << "now testing find_cn_neighbors" << std::endl;
// // std::vector<std::vector<int>> cn_neighbor = find_cn_neighbors(matrix);
// // display(cn_neighbor);
// // std::cout << "-------------------------" << std::endl;
// // std::cout << "now testing find_vn_neighbors" << std::endl;
// // std::vector<std::vector<int>> vn_neighbor = find_vn_neighbers(matrix);
// // display(vn_neighbor);