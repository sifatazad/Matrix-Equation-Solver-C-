#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define Abs(a) ((a)>=0?(a):-(a))

struct iteration {
    int *old_index;
    double *s;
    double *r;
};

static int total_variable(FILE *fp) {
    int n=0;
    char line[10000];
    while (fgets(line,sizeof(line),fp)!=NULL)
        n++;
    rewind(fp);
    return n;
}

static int compare_value(void const *a, void const *b) {
    struct iteration *it1= (struct iteration*) a;
    struct iteration *it2= (struct iteration*) b;

    if (*(it1->r) < *(it2->r))
        return 1;
    if (*(it1->r) > *(it2->r))
        return -1;
    return 0;
}

static void free_matrix(double **mat,int n) {
    for (int i=0; i<n; i++) {
        free(mat[i]);
    }
    free(mat);
}

static void solve_equation(FILE *fp, FILE *fo,double *soln) {
    int n=0;
    char line[10000];
    while (fgets(line,sizeof(line),fp)!=NULL)
        n++;
    rewind(fp);

    double **matrix=(double**)malloc(sizeof(double*)*n);
    if (!matrix) {
        return;
    }
    for (int i=0; i<n; i++) {
        matrix[i]=(double*)malloc(sizeof(double)*(n+1));
        if (!matrix[i]) {
            free_matrix(matrix,n);
            return;
        }
    }

    for (int i=0; i<n; i++) {
        for (int j=0; j<n+1; j++) {
            fscanf(fp,"%lf",&matrix[i][j]);
        }
    }

    double *scaled_vec= malloc(sizeof(double)*n);
    if (!scaled_vec) {
        free_matrix(matrix,n);
        return;
    }
    for (int i=0; i<n; i++) {
        double mxm=0;
        for (int j=0; j<n; j++) {
            if (mxm<Abs(matrix[i][j]))
                mxm=Abs(matrix[i][j]);
        }
        scaled_vec[i]=mxm;
    }

    for (int i=0; i<n; i++) {
        int it=0;
        double *ratio=malloc(sizeof(double)*(n-i));
        if (!ratio) {
            free_matrix(matrix,n);
            free(scaled_vec);
            return;
        }
        double *new_scaled=malloc(sizeof(double)*(n-i));
        if (!new_scaled) {
            free_matrix(matrix,n);
            free(scaled_vec);
            free(ratio);
            return;
        }

        for (int j=i; j<n; j++) {
            new_scaled[it++]=scaled_vec[j];
        }
        it=0;
        for (int j=i; j<n; j++) {
            ratio[it]=Abs(matrix[j][i]/new_scaled[it]);
            it++;
        }

        int *old=malloc(sizeof(int)*(n-i));
        if (!old) {
            free_matrix(matrix,n);
            free(scaled_vec);
            free(ratio);
            free(new_scaled);
            return;
        }
        for (int k=i; k<n; k++) {
            old[k-i]=k;
        }

        struct iteration *iter=malloc(sizeof(struct iteration)*(n-i));
        if (!iter) {
            free_matrix(matrix,n);
            free(scaled_vec);
            free(ratio);
            free(new_scaled);
            free(old);
            return;
        }

        for (int j=0; j<n-i; j++) {
            iter[j].r=&ratio[j];
            iter[j].s=&new_scaled[j];
            iter[j].old_index=&old[j];
        }

        qsort(iter,n-i,sizeof(struct iteration),compare_value);
        int new=*iter[0].old_index;
        double *temp=malloc(sizeof(double)*(n+1));
        if (!temp) {
            free_matrix(matrix,n);
            free(scaled_vec);
            free(ratio);
            free(new_scaled);
            free(old);
            free(iter);
            return;
        }
        for (int x=0; x<n+1; x++) {
            temp[x]=matrix[i][x];
        }
        for (int x=0; x<n+1; x++) {
            matrix[i][x]=matrix[new][x];
        }
        for (int x=0; x<n+1; x++) {
            matrix[new][x]=temp[x];
        }
        double val=matrix[i][i];
        if (Abs(val)<1e-9) {
            printf("Singular Matrix");
            free_matrix(matrix,n);
            free(scaled_vec);
            free(ratio);
            free(new_scaled);
            free(old);
            free(iter);
            free(temp);
            return;
        }
        for (int r=0; r<n+1; r++) {
            matrix[i][r]/=val;
        }

        double tmp=scaled_vec[new];
        scaled_vec[new]=scaled_vec[i];
        scaled_vec[i]=tmp;

        for (int k=0; k<n; k++) {
            if (i==k)
                continue;
            double factor= -matrix[k][i]/matrix[i][i];
            for (int col=0; col<n+1; col++) {
                matrix[k][col]+=factor*matrix[i][col];
            }
        }

        free(ratio);
        free(iter);
        free(old);
        free(temp);
        free(new_scaled);
    }
    for (int i=0; i<n; i++) {
        for (int j=0; j<n+1; j++) {
            fprintf(fo,"%.3lf\t",matrix[i][j]);
        }
        fprintf(fo,"\n");
    }
    for (int i=0; i<n; i++) {
        soln[i]=matrix[i][n];
    }
    free_matrix(matrix,n);
    free(scaled_vec);
}

static void solve_equation_without_file(const int n,double **matrix,double *soln) {
    double *scaled_vec= malloc(sizeof(double)*n);
    if (!scaled_vec) {
        return;
    }
    for (int i=0; i<n; i++) {
        double mxm=0;
        for (int j=0; j<n; j++) {
            if (mxm<Abs(matrix[i][j]))
                mxm=Abs(matrix[i][j]);
        }
        scaled_vec[i]=mxm;
    }

    for (int i=0; i<n; i++) {
        int it=0;
        double *ratio=malloc(sizeof(double)*(n-i));
        if (!ratio) {
            free(scaled_vec);
            return;
        }
        double *new_scaled=malloc(sizeof(double)*(n-i));
        if (!new_scaled) {
            free(scaled_vec);
            free(ratio);
            return;
        }

        for (int j=i; j<n; j++) {
            new_scaled[it++]=scaled_vec[j];
        }
        it=0;
        for (int j=i; j<n; j++) {
            ratio[it]=Abs(matrix[j][i]/new_scaled[it]);
            it++;
        }

        int *old=malloc(sizeof(int)*(n-i));
        if (!old) {
            free(scaled_vec);
            free(ratio);
            free(new_scaled);
            return;
        }
        for (int k=i; k<n; k++) {
            old[k-i]=k;
        }

        struct iteration *iter=malloc(sizeof(struct iteration)*(n-i));
        if (!iter) {
            free(scaled_vec);
            free(ratio);
            free(new_scaled);
            free(old);
            return;
        }

        for (int j=0; j<n-i; j++) {
            iter[j].r=&ratio[j];
            iter[j].s=&new_scaled[j];
            iter[j].old_index=&old[j];
        }

        qsort(iter,n-i,sizeof(struct iteration),compare_value);

        int new=*iter[0].old_index;
        double *temp=malloc(sizeof(double)*(n+1));
        if (!temp) {
            free(scaled_vec);
            free(ratio);
            free(new_scaled);
            free(old);
            free(iter);
            return;
        }
        for (int x=0; x<n+1; x++) {
            temp[x]=matrix[i][x];
        }
        for (int x=0; x<n+1; x++) {
            matrix[i][x]=matrix[new][x];
        }
        for (int x=0; x<n+1; x++) {
            matrix[new][x]=temp[x];
        }
        double val=matrix[i][i];
        if (Abs(val)<1e-9) {
            printf("Singular Matrix");
            free(scaled_vec);
            free(ratio);
            free(new_scaled);
            free(old);
            free(iter);
            free(temp);
            return;
        }
        for (int r=0; r<n+1; r++) {
            matrix[i][r]/=val;
        }

        double tmp=scaled_vec[new];
        scaled_vec[new]=scaled_vec[i];
        scaled_vec[i]=tmp;

        for (int k=0; k<n; k++) {
            if (i==k)
                continue;
            double factor= -matrix[k][i]/matrix[i][i];
            for (int col=0; col<n+1; col++) {
                matrix[k][col]+=factor*matrix[i][col];
            }
        }

        free(ratio);
        free(iter);
        free(old);
        free(temp);
        free(new_scaled);
    }
    for (int i=0; i<n; i++) {
        soln[i]=matrix[i][n];
    }
    // free_matrix(matrix,n);
    free(scaled_vec);
}

static void inverse_matrix(FILE *fp,FILE *fo, double **ans) {
    int n=0;
    double det=1;
    char line[10000];
    while (fgets(line,sizeof(line),fp)!=NULL)
        n++;
    rewind(fp);

    double **matrix=(double**)malloc(sizeof(double*)*n);
    if (!matrix) {
        return;
    }
    for (int i=0; i<n; i++) {
        matrix[i]=(double*)malloc(sizeof(double)*n);
        if (!matrix[i]) {
            free_matrix(matrix,n);
            return;
        }
    }

    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            fscanf(fp,"%lf",&matrix[i][j]);
        }
    }

    double **eye=(double**)malloc(sizeof(double*)*n);
    if (!eye) {
        free_matrix(matrix,n);
        return;
    }
    for (int i=0; i<n; i++) {
        eye[i]=(double*)malloc(sizeof(double)*n);
        if (!eye[i]) {
            free_matrix(matrix,n);
            free_matrix(eye,n);
            return;
        }
    }

    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            if (i==j)
                eye[i][j]=1;
            else
                eye[i][j]=0;
        }
    }

    double **aug=(double**)malloc(sizeof(double*)*n);
    if (!aug) {
        free_matrix(matrix,n);
        free_matrix(eye,n);
        return;
    }
    for (int i=0; i<n; i++) {
        aug[i]=(double*)malloc(sizeof(double)*(n+1));
        if (!aug[i]) {
            free_matrix(matrix,n);
            free_matrix(eye,n);
            free_matrix(aug,n);
            return;
        }
    }

    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            for (int k=0; k<n+1; k++) {
                if (k!=n) {
                    aug[j][k]=matrix[j][k];
                }
                else {
                    aug[j][k]=eye[j][i];
                }
            }
        }
        double *small_ans= malloc(sizeof(double)*n);
        if (!small_ans) {
            free_matrix(matrix,n);
            free_matrix(eye,n);
            free_matrix(aug,n);
            return;
        }
        solve_equation_without_file(n,aug,small_ans);
        int it=0;
        for (int j=0; j<n; j++) {
            ans[j][i]=small_ans[it++];
        }
        free(small_ans);
    }
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            fprintf(fo,"%.3lf\t",ans[i][j]);
        }
        fprintf(fo,"\n");
    }
    free_matrix(matrix,n);
    free_matrix(aug,n);
    free_matrix(eye,n);
}

static void inverse_matrix_without_file(const int n,double **matrix,double **ans) {
    double **eye=(double**)malloc(sizeof(double*)*n);
    if (!eye) {
        return;
    }
    for (int i=0; i<n; i++) {
        eye[i]=(double*)malloc(sizeof(double)*n);
        if (!eye[i]) {
            free_matrix(eye,n);
            return;
        }
    }

    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            if (i==j)
                eye[i][j]=1;
            else
                eye[i][j]=0;
        }
    }

    double **aug=(double**)malloc(sizeof(double*)*n);
    if (!aug) {
        free_matrix(eye,n);
        return;
    }
    for (int i=0; i<n; i++) {
        aug[i]=(double*)malloc(sizeof(double)*(n+1));
        if (!aug[i]) {
            free_matrix(eye,n);
            free_matrix(aug,n);
            return;
        }
    }

    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            for (int k=0; k<n+1; k++) {
                if (k!=n) {
                    aug[j][k]=matrix[j][k];
                }
                else {
                    aug[j][k]=eye[j][i];
                }
            }
        }
        double *small_ans= malloc(sizeof(double)*n);
        if (!small_ans) {
            free_matrix(eye,n);
            free_matrix(aug,n);
            return;
        }
        solve_equation_without_file(n,aug,small_ans);
        int it=0;
        for (int j=0; j<n; j++) {
            ans[j][i]=small_ans[it++];
        }
        free(small_ans);
    }
    free_matrix(aug,n);
    free_matrix(eye,n);
}

static double determinant(FILE *fp) {
    int n=0;
    double det=1;
    char line[10000];
    while (fgets(line,sizeof(line),fp)!=NULL)
        n++;
    rewind(fp);

    double **matrix=(double**)malloc(sizeof(double*)*n);
    if (!matrix) {
        return NAN;
    }
    for (int i=0; i<n; i++) {
        matrix[i]=(double*)malloc(sizeof(double)*n);
        if (!matrix[i]) {
            free_matrix(matrix,n);
            return NAN;
        }
    }

    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            fscanf(fp,"%lf",&matrix[i][j]);
        }
    }

    double *scaled_vec= malloc(sizeof(double)*n);
    if (!scaled_vec) {
        free_matrix(matrix,n);
        return NAN;
    }
    for (int i=0; i<n; i++) {
        double mxm=0;
        for (int j=0; j<n; j++) {
            if (mxm<Abs(matrix[i][j]))
                mxm=Abs(matrix[i][j]);
        }
        scaled_vec[i]=mxm;
    }

    for (int i=0; i<n; i++) {
        int it=0;
        double *ratio=malloc(sizeof(double)*(n-i));
        if (!ratio) {
            free_matrix(matrix,n);
            free(scaled_vec);
            return NAN;
        }
        double *new_scaled=malloc(sizeof(double)*(n-i));
        if (!new_scaled) {
            free_matrix(matrix,n);
            free(scaled_vec);
            free(ratio);
            return NAN;
        }

        for (int j=i; j<n; j++) {
            new_scaled[it++]=scaled_vec[j];
        }
        it=0;
        for (int j=i; j<n; j++) {
            ratio[it]=Abs(matrix[j][i]/new_scaled[it]);
            it++;
        }

        int *old=malloc(sizeof(int)*(n-i));
        if (!old) {
            free_matrix(matrix,n);
            free(scaled_vec);
            free(ratio);
            free(new_scaled);
            return NAN;
        }
        for (int k=i; k<n; k++) {
            old[k-i]=k;
        }

        struct iteration *iter=malloc(sizeof(struct iteration)*(n-i));
        if (!iter) {
            free_matrix(matrix,n);
            free(scaled_vec);
            free(ratio);
            free(new_scaled);
            free(old);
            return NAN;
        }

        for (int j=0; j<n-i; j++) {
            iter[j].r=&ratio[j];
            iter[j].s=&new_scaled[j];
            iter[j].old_index=&old[j];
        }

        qsort(iter,n-i,sizeof(struct iteration),compare_value);
        int new=*iter[0].old_index;
        double *temp=malloc(sizeof(double)*n);
        if (!temp) {
            free_matrix(matrix,n);
            free(scaled_vec);
            free(ratio);
            free(new_scaled);
            free(old);
            free(iter);
            return NAN;
        }
        for (int x=0; x<n; x++) {
            temp[x]=matrix[i][x];
        }
        for (int x=0; x<n; x++) {
            matrix[i][x]=matrix[new][x];
        }
        for (int x=0; x<n; x++) {
            matrix[new][x]=temp[x];
        }
        double val=matrix[i][i];
        det*=val;
        if (new!=i) {
            det*=-1;
        }
        if (Abs(val)<1e-9) {
            printf("Singular Matrix");
            free_matrix(matrix,n);
            free(scaled_vec);
            free(ratio);
            free(new_scaled);
            free(old);
            free(iter);
            free(temp);
            return 0;
        }
        for (int r=0; r<n; r++) {
            matrix[i][r]/=val;
        }

        double tmp=scaled_vec[new];
        scaled_vec[new]=scaled_vec[i];
        scaled_vec[i]=tmp;

        for (int k=0; k<n; k++) {
            if (i==k)
                continue;
            double factor= -matrix[k][i]/matrix[i][i];
            for (int col=0; col<n; col++) {
                matrix[k][col]+=factor*matrix[i][col];
            }
        }

        free(ratio);
        free(iter);
        free(old);
        free(temp);
        free(new_scaled);
    }

    free_matrix(matrix,n);
    free(scaled_vec);
    return det;
}

static double determinant_without_file(const int n,double **matrix) {
    double det=1;

    double *scaled_vec= malloc(sizeof(double)*n);
    if (!scaled_vec) {
        free_matrix(matrix,n);
        return NAN;
    }
    for (int i=0; i<n; i++) {
        double mxm=0;
        for (int j=0; j<n; j++) {
            if (mxm<Abs(matrix[i][j]))
                mxm=Abs(matrix[i][j]);
        }
        scaled_vec[i]=mxm;
    }

    for (int i=0; i<n; i++) {
        int it=0;
        double *ratio=malloc(sizeof(double)*(n-i));
        if (!ratio) {
            free_matrix(matrix,n);
            free(scaled_vec);
            return NAN;
        }
        double *new_scaled=malloc(sizeof(double)*(n-i));
        if (!new_scaled) {
            free_matrix(matrix,n);
            free(scaled_vec);
            free(ratio);
            return NAN;
        }

        for (int j=i; j<n; j++) {
            new_scaled[it++]=scaled_vec[j];
        }
        it=0;
        for (int j=i; j<n; j++) {
            ratio[it]=Abs(matrix[j][i]/new_scaled[it]);
            it++;
        }

        int *old=malloc(sizeof(int)*(n-i));
        if (!old) {
            free_matrix(matrix,n);
            free(scaled_vec);
            free(ratio);
            free(new_scaled);
            return NAN;
        }
        for (int k=i; k<n; k++) {
            old[k-i]=k;
        }

        struct iteration *iter=malloc(sizeof(struct iteration)*(n-i));
        if (!iter) {
            free_matrix(matrix,n);
            free(scaled_vec);
            free(ratio);
            free(new_scaled);
            free(old);
            return NAN;
        }

        for (int j=0; j<n-i; j++) {
            iter[j].r=&ratio[j];
            iter[j].s=&new_scaled[j];
            iter[j].old_index=&old[j];
        }

        qsort(iter,n-i,sizeof(struct iteration),compare_value);
        int new=*iter[0].old_index;
        double *temp=malloc(sizeof(double)*n);
        if (!temp) {
            free_matrix(matrix,n);
            free(scaled_vec);
            free(ratio);
            free(new_scaled);
            free(old);
            free(iter);
            return NAN;
        }
        for (int x=0; x<n; x++) {
            temp[x]=matrix[i][x];
        }
        for (int x=0; x<n; x++) {
            matrix[i][x]=matrix[new][x];
        }
        for (int x=0; x<n; x++) {
            matrix[new][x]=temp[x];
        }
        double val=matrix[i][i];
        det*=val;
        if (new!=i) {
            det*=-1;
        }
        if (Abs(val)<1e-9) {
            printf("Singular Matrix");
            free_matrix(matrix,n);
            free(scaled_vec);
            free(ratio);
            free(new_scaled);
            free(old);
            free(iter);
            free(temp);
            return 0;
        }
        for (int r=0; r<n; r++) {
            matrix[i][r]/=val;
        }

        double tmp=scaled_vec[new];
        scaled_vec[new]=scaled_vec[i];
        scaled_vec[i]=tmp;

        for (int k=0; k<n; k++) {
            if (i==k)
                continue;
            double factor= -matrix[k][i]/matrix[i][i];
            for (int col=0; col<n; col++) {
                matrix[k][col]+=factor*matrix[i][col];
            }
        }

        free(ratio);
        free(iter);
        free(old);
        free(temp);
        free(new_scaled);
    }

    free_matrix(matrix,n);
    free(scaled_vec);
    return det;
}

int main() {
    FILE *fp=fopen("matrix.txt","r");
    if (!fp)
        return 1;
    FILE *fo= fopen("output.txt","w");
    if (!fo) {
        fclose(fp);
        return 1;
    }

    int n=total_variable(fp);
    double *solution= malloc(sizeof(double)*n);
    if (!solution) {
        fclose(fp);
        fclose(fo);
        return 1;
    }
    solve_equation(fp,fo,solution);
    // for (int i=0; i<n; i++) {
    //     printf("%.3lf\n",solution[i]);
    // }

    FILE *fp1=fopen("sample_matrix.txt","r");
    if (!fp1) {
        fclose(fp);
        fclose(fo);
        free(solution);
        return 1;
    }
    n=total_variable(fp1);
    double **ans=malloc(sizeof(double*)*n);
    if (!ans) {
        fclose(fp);
        fclose(fo);
        free(solution);
        return 1;
    }
    for (int i=0; i<n; i++) {
        ans[i]=malloc(sizeof(double)*n);
        if (!ans[i]) {
            fclose(fp);
            fclose(fo);
            free(solution);
            free_matrix(ans,n);
            return 1;
        }
    }
    FILE *fo1=fopen("inv_matrix.txt","w");
    if (!fo1) {
        fclose(fp);
        fclose(fo);
        free(solution);
        free_matrix(ans,n);;
        return 1;
    }
    inverse_matrix(fp1,fo1,ans);
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            printf("%.3lf\t",ans[i][j]);
        }
        printf("\n");
    }

    FILE *fp2=fopen("det.txt","r");
    if (!fp2) {
        fclose(fp);
        fclose(fo);
        free(solution);
        free_matrix(ans,n);
        fclose(fo1);
        return 1;
    }
    printf("%.3lf\n",determinant(fp2));

    fclose(fp);
    fclose(fo);
    fclose(fp1);
    fclose(fo1);
    fclose(fp2);
    free(solution);
    free_matrix(ans,n);
    return 0;
}