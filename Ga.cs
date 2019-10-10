using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeneticAlgorithm
{
    class Ga
    {
        private int scale;                    // 种群规模
        private int cityNum;                  // 城市数量，染色体长度
        private int MAX_GEN;                  // 运行代数
        private double[,] distance;           // 距离矩阵
        private int bestGen;                  // 最佳出现代数
        private double bestDistance;          // 最佳长度
        private int[] bestRoute;              // 最佳路径

 
        private int[,] oldPopulation;
        private int[,] newPopulation;         // 新的种群
        private double[] fitness;             // 种群适应度

        private float[] Pi;                   // 种群中各个个体的累计概率
        private float Pc;                     // 交叉概率
        private float Pm;                     // 变异概率
        private int Gen;                      // 当前代数

        private int V = 10;                   //速度
        private double[] idealDistance;       //理想距离

        private Random random;

        public Ga()
        {

        }


        public Ga(int s, int n, int g, float c, float m)
        {
            scale = s;
            cityNum = n;
            MAX_GEN = g;
            Pc = c;
            Pm = m;
        }

        /// <summary>
        /// 初始化
        /// </summary>
        /// <param name="filename">数据存放的地址</param>
        public void Init(string filename)
        {
            int[] x;
            int[] y;
            string strbuff;

            using (StreamReader data = new StreamReader(filename, Encoding.Default))
            {

                distance = new double[cityNum, cityNum];

                x = new int[cityNum];
                y = new int[cityNum];

                for (int i = 0; i < cityNum; i++)
                {
                    strbuff = data.ReadLine();

                    // 字符分割
                    string[] strcol = strbuff.Split();

                    x[i] = Convert.ToInt32(strcol[1]);  // x坐标
                    y[i] = Convert.ToInt32(strcol[2]);  // y坐标
                }
            }

            // 计算距离矩阵			
            for (int i = 0; i < cityNum - 1; i++)
            {
                distance[i, i] = 0;  
                for (int j = i + 1; j < cityNum; j++)
                {
                    double rij = Math
                            .Sqrt(((x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j])
                                    * (y[i] - y[j])) / 10.0);
                    // 四舍五入，取整
                    int tij = (int)Math.Round(rij);
                    if (tij < rij)
                    {
                        distance[i, j] = tij + 1;
                        distance[j, i] = distance[i, j];
                    }
                    else
                    {
                        distance[i, j] = tij;
                        distance[j, i] = distance[i, j];
                    }
                }
            }
            distance[cityNum - 1, cityNum - 1] = 0;

            bestDistance = int.MaxValue;
            bestRoute = new int[cityNum + 1];
            bestGen = 0;
            Gen = 0;

            newPopulation = new int[scale, cityNum];
            oldPopulation = new int[scale, cityNum];
            fitness = new double[scale];
            Pi = new float[scale];

            random = new Random(DateTime.Now.Millisecond);

            //随机生成限制条件
            Random randomNum = new Random();
            idealDistance = new double[cityNum];
            for (int r = 0; r < cityNum; r++)
            {
                if (r == 10 || r == 20 || r == 30 || r== 38 )
                {
                    idealDistance[r] = r * randomNum.Next(40,60) * V;
                }
                else
                {
                    idealDistance[r] = int.MaxValue;
                }
                
            }

        }

        /// <summary>
        /// 初始化种群
        /// </summary>
        private void InitGroup()
        {
            int i, j, k;
            for (k = 0; k < scale; k++)
            {
                oldPopulation[k, 0] = random.Next(65535) % cityNum;
                for (i = 1; i < cityNum;)
                {
                    oldPopulation[k, i] = random.Next(65535) % cityNum;
                    for (j = 0; j < i; j++)
                    {
                        if (oldPopulation[k, i] == oldPopulation[k, j])
                        {
                            break;
                        }
                    }
                    if (j == i)
                    {
                        i++;
                    }

                }
            }
        }

        /// <summary>
        /// 计算适应度函数值
        /// </summary>
        /// <param name="chromosome"></param>
        /// <returns></returns>
        private double Evaluate(int[] chromosome)
        {
            double len = 0.0;
            double countDiff = 0.0;    //累计差值
            double penatly = 0.0;      //惩罚函数

            // 染色体，起始城市,城市1,城市2...城市n
            for (int i = 1; i < cityNum; i++)
            {
                len += distance[chromosome[i - 1], chromosome[i]];
                
                Getpenatly(len,ref countDiff, ref penatly,chromosome[i]);
            }

            // 城市n,起始城市
            len += (distance[chromosome[cityNum - 1], chromosome[0]] + penatly);

            return len;
        }


        /// <summary>
        /// 遍历每个节点，是否超过理想的距离
        /// </summary>
        private void Getpenatly(double len,ref double countDiff, ref double penatly,int cityNode)
        {         
            double idealDis = idealDistance[cityNode];
            if (len > idealDis)
            {
                countDiff += len - idealDis;
                penatly += Math.Pow(0.1 * Gen, 0.8) * Math.Pow(countDiff, 0.8);
              
            }//大于理想距离
            
        }
      
        /// <summary>
        /// 计算种群中各个个体的累积概率
        /// </summary>
        private void CountRate()
        {
            int k;
            double sumFitness = 0;

            double[] tempf = new double[scale];

            for (k = 0; k < scale; k++)
            {
                tempf[k] = 10.0 / fitness[k];
                sumFitness += tempf[k];
            }

            Pi[0] = (float)(tempf[0] / sumFitness);
            for (k = 1; k < scale; k++)
            {
                Pi[k] = (float)(tempf[k] / sumFitness + Pi[k - 1]);
            }

        }

        
        /// <summary>
        /// 挑选某代种群中适应度最高的个体，直接复制到子代中
        /// </summary>
        private void SelectBestIndividual()
        {
            int k, i, maxId;
            double maxEvaluation;

            maxId = 0;
            maxEvaluation = fitness[0];
            for (k = 1; k < scale; k++)
            {
                if (maxEvaluation > fitness[k])
                {
                    maxEvaluation = fitness[k];
                    maxId = k;
                }
            }

            if (bestDistance > maxEvaluation)
            {
                bestDistance = maxEvaluation;
                bestGen = Gen;
                for (i = 0; i < cityNum; i++)
                {
                    bestRoute[i] = oldPopulation[maxId, i];
                }
            }

            ProductIndividual(0, maxId);// 将当代种群中适应度最高的染色体k复制到新种群中，排在第一位0
        }

        /// <summary>
        /// 复制染色体
        /// </summary>
        /// <param name="k">新染色体在种群中的位置</param>
        /// <param name="kk">旧的染色体在种群中的位置</param>
        private void ProductIndividual(int k, int kk)
        {
            int i;
            for (i = 0; i < cityNum; i++)
            {
                newPopulation[k, i] = oldPopulation[kk, i];
            }
        }

        /// <summary>
        /// 赌轮选择策略挑选
        /// </summary>
        private void Select()
        {
            int k, i, selectId;
            float ran1;
            for (k = 1; k < scale; k++)
            {
                ran1 = (float)(random.Next(65535) % 1000 / 1000.0);
                // 产生方式
                for (i = 0; i < scale; i++)
                {
                    if (ran1 <= Pi[i])
                    {
                        break;
                    }
                }
                selectId = i;
 
                ProductIndividual(k, selectId);
            }
        }

        /// <summary>
        /// 进化函数，正常交叉变异
        /// </summary>
        private void Evolution()
        {
            int k;
            // 挑选某代种群中适应度最高的个体
            SelectBestIndividual();

            // 赌轮选择策略挑选scale-1个下一代个体
            Select();

            double r;

            // 交叉及变异
            for (k = 0; k < scale; k = k + 2)
            {
                r = random.NextDouble();
                if (r < Pc)
                {
                    OXCross1(k, k + 1);
                }
                else
                {
                    r = random.NextDouble();
                    if (r < Pm)
                    {
                        OnCVariation(k);
                    }
                    r = random.NextDouble();
                    if (r < Pm)
                    {
                        OnCVariation(k + 1);
                    }
                }

            }
        }

        /// <summary>
        /// 进化函数，保留最好染色体不进行交叉变异
        /// </summary>
        private void Evolution1()
        {
            int k;
            // 挑选某代种群中适应度最高的个体
            SelectBestIndividual();

            // 赌轮选择策略挑选scale-1个下一代个体
            Select();

            double r;

            for (k = 1; k + 1 < scale / 2; k = k + 2)
            {
                r = random.NextDouble();
                if (r < Pc)
                {
                    OXCross1(k, k + 1);
                }
                else
                {
                    r = random.NextDouble();
                    if (r < Pm)
                    {
                        OnCVariation(k);
                    }
                    r = random.NextDouble();
                    if (r < Pm)
                    {
                        OnCVariation(k + 1);
                    }
                }
            }
            if (k == scale / 2 - 1)// 剩最后一个染色体没有交叉L-1
            {
                r = random.NextDouble();
                if (r < Pm)
                {
                    OnCVariation(k);
                }
            }

        }

        /// <summary>
        /// 顺序交叉
        /// </summary>
        /// <param name="k1">交叉的染色体序号</param>
        /// <param name="k2">交叉的染色体序号</param>
        private void OXCross(int k1, int k2)
        {
            int i, j, k, flag;
            int ran1, ran2, temp;
            int[] Gh1 = new int[cityNum];
            int[] Gh2 = new int[cityNum];


            ran1 = random.Next(65535) % cityNum;   //交叉的起始位置
            ran2 = random.Next(65535) % cityNum;   //交叉的结束位置

            while (ran1 == ran2)
            {
                ran2 = random.Next(65535) % cityNum;
            }

            if (ran1 > ran2)// 确保ran1<ran2
            {
                temp = ran1;
                ran1 = ran2;
                ran2 = temp;
            }


            flag = ran2 - ran1 + 1;// 删除重复基因前染色体长度
            for (i = 0, j = ran1; i < flag; i++, j++)
            {
                Gh1[i] = newPopulation[k2, j];
                Gh2[i] = newPopulation[k1, j];
            }
            // 已近赋值i=ran2-ran1个基因

            for (k = 0, j = flag; j < cityNum;)

            {
                Gh1[j] = newPopulation[k1, k++];
                for (i = 0; i < flag; i++)
                {
                    if (Gh1[i] == Gh1[j])
                    {
                        break;
                    }
                }
                if (i == flag)
                {
                    j++;
                }
            }

            for (k = 0, j = flag; j < cityNum;)
            {
                Gh2[j] = newPopulation[k2, k++];
                for (i = 0; i < flag; i++)
                {
                    if (Gh2[i] == Gh2[j])
                    {
                        break;
                    }
                }
                if (i == flag)
                {
                    j++;
                }
            }

            for (i = 0; i < cityNum; i++)
            {
                newPopulation[k1, i] = Gh1[i];
                newPopulation[k2, i] = Gh2[i];// 交叉完毕放回种群
            }

        }

        /// <summary>
        /// 交叉操作,相同染色体交叉产生不同子代染色体
        /// </summary>
        /// <param name="k1"></param>
        /// <param name="k2"></param>
        private void OXCross1(int k1, int k2)
        {
            int i, j, k, flag;
            int ran1, ran2, temp;
            int[] Gh1 = new int[cityNum];
            int[] Gh2 = new int[cityNum];

            ran1 = random.Next(65535) % cityNum;
            ran2 = random.Next(65535) % cityNum;
            while (ran1 == ran2)
            {
                ran2 = random.Next(65535) % cityNum;
            }

            if (ran1 > ran2)// 确保ran1<ran2
            {
                temp = ran1;
                ran1 = ran2;
                ran2 = temp;
            }

            // 将染色体1中的第三部分移到染色体2的首部
            for (i = 0, j = ran2; j < cityNum; i++, j++)
            {
                Gh2[i] = newPopulation[k1, j];
            }

            flag = i;// 染色体2原基因开始位置

            for (k = 0, j = flag; j < cityNum;)
            {
                Gh2[j] = newPopulation[k2, k++];
                for (i = 0; i < flag; i++)
                {
                    if (Gh2[i] == Gh2[j])
                    {
                        break;
                    }
                }
                if (i == flag)
                {
                    j++;
                }
            }

            flag = ran1;
            for (k = 0, j = 0; k < cityNum;)
            {
                Gh1[j] = newPopulation[k1, k++];
                for (i = 0; i < flag; i++)
                {
                    if (newPopulation[k2, i] == Gh1[j])
                    {
                        break;
                    }
                }
                if (i == flag)
                {
                    j++;
                }
            }

            flag = cityNum - ran1;

            for (i = 0, j = flag; j < cityNum; j++, i++)
            {
                Gh1[j] = newPopulation[k2, i];
            }

            for (i = 0; i < cityNum; i++)
            {
                newPopulation[k1, i] = Gh1[i];
                newPopulation[k2, i] = Gh2[i];// 交叉完毕放回种群
            }
        }

        /// <summary>
        /// 变异 多次对换变异算子
        /// </summary>
        /// <param name="k"></param>
        private void OnCVariation(int k)
        {
            int ran1, ran2, temp;
            int count;// 对换次数

            count = random.Next(65535) % cityNum;

            for (int i = 0; i < count; i++)
            {

                ran1 = random.Next(65535) % cityNum;
                ran2 = random.Next(65535) % cityNum;
                while (ran1 == ran2)
                {
                    ran2 = random.Next(65535) % cityNum;
                }
                temp = newPopulation[k, ran1];
                newPopulation[k, ran1] = newPopulation[k, ran2];
                newPopulation[k, ran2] = temp;
            }


        }
   
        /// <summary>
        /// 开始
        /// </summary>
        public void Solve()
        {
            int i;
            int k;

            // 初始化种群
            InitGroup();
            // 计算初始化种群适应度
            for (k = 0; k < scale; k++)
            {
                int[] chrom = new int[cityNum];
                for (int c = 0; c < cityNum; c++)
                {
                    chrom[c] = oldPopulation[k, c];
                }
                fitness[k] = Evaluate(chrom);
                
            }
            CountRate();
            
            for (Gen = 0; Gen < MAX_GEN; Gen++)
            {
                //evolution1();
                Evolution();
                // 将新种群newGroup复制到旧种群oldGroup中，准备下一代进化
                for (k = 0; k < scale; k++)
                {
                    for (i = 0; i < cityNum; i++)
                    {
                        oldPopulation[k, i] = newPopulation[k, i];
                    }
                }
                // 计算种群适应度
                for (k = 0; k < scale; k++)
                {
                    int[] chrom = new int[cityNum];
                    for (int c = 0; c < cityNum; c++)
                    {
                        chrom[c] = oldPopulation[k, c];
                    }
                    fitness[k] = Evaluate(chrom);
                }
                // 计算种群中各个个体的累积概率
                CountRate();
            }

            Console.WriteLine("最后种群...");
            for (k = 0; k < scale; k++)
            {
                for (i = 0; i < cityNum; i++)
                {
                    Console.Write(oldPopulation[k, i] + ",");
                }
                Console.WriteLine();
                Console.WriteLine("---" + fitness[k] + " " + Pi[k]);
            }

            Console.WriteLine("最佳长度出现代数：");
            Console.WriteLine(bestGen);
            Console.WriteLine("最佳长度");
            Console.WriteLine(bestDistance);
            Console.WriteLine("最佳路径：");

            for (i = 0; i < cityNum; i++)
            {
                Console.Write(bestRoute[i] + ",");
            }

        }
    }
}


