typedef struct
{
	int N; /*隐藏状态数目; Q={1,2,……N}*/
	int M; /*观测符号数目; V={1,2,……M}*/
	double **A; /*状态转移矩阵，是一个二维数组 
	A[1……N][1……N],a[i][j]是从t时刻状态i到t+1时刻状态j的转移概率*/
	double **B; /*混淆矩阵(观测概率矩阵)B[1……N][1……M], b[j][k]在隐藏状态j时生成观测符号K的概率*/
	double *pi; /*初始向量pi[1……N], pi是初始状态概率分布*/
}HMM;


/**
函数参数说明：
*phmm:已知的HMM模型，T:观测符号序列长度
*O:观测序列；**alpha：局部概率；*pprob:最终的观测概率
*/
void forward(HMM *phmm, int T, int *O, double **alpha, double **pprob)
{
	int i,j; 	/*状态索引*/
	int t;	/*时间索引*/
	double sum; /*求局部概率时的中间值*/
	
	/*1.初始化，计算t=1时刻所有状态的局部概率alpha:
		t=1时刻，观测符号是确定的
		就是状态i自己的概率 * 混淆矩阵中，i生成观测符号的概率
	*/
	for(i=1, i<=phmm->N; i++)
		alpha[1][i]=phmm->pi[i]*phmm->B[i][O[1]];  
	
	/*2.归纳：递归计算每个时间点，t=2,……T时候的局部概率
	就是那个求和公式，
	*/
	for(t=1; t<T; t++) /*每一个时刻t*/
	{
		for(j=1; j<=phmm->N; j++) /*t+1时刻任何一个隐藏状态j的局部概率的计算*/
		{
			sum=0.0;/*t+1时刻任何一个隐藏状态j生成观测符号的局部概率都需要前面t时刻所有的隐层状态做基础*/
			for(i=1; i<=phmm->N; i++) /*t时刻每个状态i本身的概率 乘以 状态转移矩阵*/
				sum += alpha[t][i]*(phmm->A[i][j])
			alpha[t+1][j]=sum*(phmm->B[j][O[T+1]])
		}
	}
	/* 3. 终止：观察序列的概率等于T时刻所有局部概率之和*/
	*pprob=0.0;
	for(i=1; i<phmm->N; i++)
		*pprob +=alpha[T][i];
	
}

void Viterbi(HMM *phmm, int T, int *O, double **delta, int **psi, int *q, double *pprob)
{
	int i, j; /*state indices */
	int t;	/*time index */
	int maxvalind;
	double maxval, val;
	
	/*初始化 */
	for(i=1; i<=phmm->N; i++)
	{
		delta[1][i]=phmm->pi[i]*(phmm->B[i][O[1]]);
		psi[1][i]=0;
	}
	
	/* recursion 递归 一个是得到局部概率中的最优的，一个是反向指针 */
	for(t=2; t<=T; t++)
	{
		for(j=1; j<=phmm->N; j++)
		{
			maxval=0.0;
			maxvalind = 1;
			for(i=1; i<=phmm->N; i++)
			{
				val = delta[t-1][i]*(phmm->A[i][j]);
				if(val>maxval)
				{
					maxval=val;
					maxvalind=i;
				}
			}
			delta[t][j]=maxval*(phmm->B[j][O[t]]);
			psi[t][j]=maxvalind;
		}
	}
	
	/* Termination */
	*pprob=0.0;
	q[T]=1;
	for(i=1; i<=phmm->N; i++)
	{
		if(delta[T][i]>*pprob)
		{
			*pprob=delta[T][i]
			q[T]=i;
		}
	}
	
	/* Path(state sequence) backtracking */
	for(t=T-1; t>=1; t--)
		q[t]=psi[T+1][q[t+1]];
	
}


