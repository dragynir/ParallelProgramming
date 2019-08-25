#include <pthread.h>
#include <stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<time.h>
#include<deque>
#include<set>
//for sleap
#include <chrono>
#include <thread>

//#define DEBUG
//#define LVL_DEBUG

using std::deque;
using std::set;
const int PROC_OUT = -1;

typedef struct Task{
	int time_mls;
}task;
//300 - 500
//Hight balance
//140 tasks
//75410

//Hight hard
//140
//562157


//Medium balance
//140
//39738

//Medium hard
//140
//64696


//Low balance
//140
//41312

//Low hard
//140
//38920




//mpicxx dynamictasks.cpp -o a.out -std=c++11 -O3 -lpthread -mt_mpi


//std::this_thread::sleep_for(std::chrono::milliseconds(x));

const int TAG_GET = 12;
const int TAG_ANS = 21;
int p_count, m_rank;
int lvl;
pthread_mutex_t mutex;
pthread_t worker, listener, observer;

set<int> proc;
set<int> unAnsDied;
deque<task> tasks;


int PROVIDED_TASKS;
const int MIN_V = 300, MAX_V = 500;

int globalWorkingTime;//будет учитываться позже во время балансировки задач
int finishedTasks;
bool workEnd = false;



/////////////////////////////////////////INIT
void initTasksLow(int iter){
	task t;
	int myTasksCount = PROVIDED_TASKS / p_count;
	if(m_rank < PROVIDED_TASKS%p_count){
		++myTasksCount;
	}
	for (int i = 0; i < myTasksCount; ++i)
	{
		t.time_mls = MIN_V + rand() % (MAX_V - MIN_V);
		tasks.push_back(t);
	}
	//printf("Rank %d , %d\n", m_rank, myTasksCount);
}



void initTasksMedium(int iter){
	if(p_count == 1){
		initTasksLow(iter);
		return;
	}
	task t;
	int myTasksCount = 0;

	int d = (2 * PROVIDED_TASKS) / (p_count * (p_count - 1));
	int lost = PROVIDED_TASKS - d * p_count * (p_count - 1) / 2;
	if(m_rank < lost){
		++myTasksCount;
	}
	if(m_rank < lost % p_count){
		++myTasksCount;
	}
	myTasksCount+=(d * m_rank);

	for (int i = 0; i < myTasksCount; ++i)
	{
		t.time_mls = MIN_V + rand() % (MAX_V - MIN_V);
		tasks.push_back(t);
	}
	//printf("Rank %d , %d\n", m_rank, myTasksCount);
}


void initTasksHigh(int iter){
	task t;

	if(iter % p_count == m_rank){
		for (int i = 0; i < PROVIDED_TASKS; ++i)
		{
			t.time_mls = MIN_V + rand() % (MAX_V - MIN_V);
			tasks.push_back(t);
		}
		//printf("Start rang: %d , %d\n", m_rank, PROVIDED_TASKS);
	}
}



void printTasks(){
	for (int i = 0; i < tasks.size(); ++i)
	{
		//printf("Task %d: %d mls \n", i, tasks[i].time_mls);
	}
}
/////////////////////////////////////////INIT





/////////////////////////////////////////Worker

void* startWork(void* t){
	//deque<task>& tasks = *((deque<task>*)t);
	while(true){
		task* cur_task = nullptr;

		pthread_mutex_lock(&mutex);
		if(workEnd){
			pthread_mutex_unlock(&mutex);
			break;
		}else if(lvl == 1 && tasks.empty()){
			pthread_mutex_unlock(&mutex);
			break;
		}  

		if(!tasks.empty()){
			//printf("%d:%d\n",m_rank, tasks.size());
			cur_task = &tasks.front();
			tasks.pop_front();
		}
		pthread_mutex_unlock(&mutex);

		if(cur_task != nullptr){
			std::this_thread::sleep_for(std::chrono::milliseconds(cur_task->time_mls));
			pthread_mutex_lock(&mutex);
			//globalWorkingTime+=cur_task->time_mls;
			++finishedTasks;
			pthread_mutex_unlock(&mutex);
		}
	}
	#ifdef DEBUG
	printf("%d , Worker end\n", m_rank);
	#endif
}




/////////////////////////////////////////Worker



/////////////////////////////////////////Observer




void* startObserv(void* t){
	pthread_mutex_lock(&mutex);
	if(lvl == 1){
		pthread_mutex_unlock(&mutex);

		#ifdef LVL_DEBUG
		printf("Observer lvl 1 end.\n");
		#endif

		return 0;
	} 
	pthread_mutex_unlock(&mutex);

	MPI_Status status;
	set<int>::iterator it;
	while(true){
		int tasks_count = 0;
		//можно убрать этот блок
		pthread_mutex_lock(&mutex);
		if(workEnd){
			pthread_mutex_unlock(&mutex);
			break;
		}
		tasks_count = tasks.size();
		pthread_mutex_unlock(&mutex);
		//можно убрать этот блок

		if(tasks_count < PROVIDED_TASKS / p_count){
			set<int> l_proc(proc);
		
			for (it=l_proc.begin(); it!=l_proc.end(); ++it){
    			int i = *it;
				int res = 0;
				//будет передавать кол-во ост задач

				MPI_Send(&tasks_count, 1, MPI_INT, i, TAG_GET, MPI_COMM_WORLD);
				MPI_Recv(&res, 1, MPI_INT, i, TAG_ANS, MPI_COMM_WORLD, &status);
				
				if(res >= 0){
					if(res == 0) continue;
					task new_task;
					new_task.time_mls = res;
					pthread_mutex_lock(&mutex);
					tasks.push_back(new_task);
					pthread_mutex_unlock(&mutex);
					break;
				}else if(res == PROC_OUT){
					#ifdef DEBUG
					printf("Recv OUT\n");
					#endif
					proc.erase(i);
				}
			}//for


			//за время опроса могли кончиться задачи
			//поэтому проверим что их нет 
			//для того чтобы не проводить еще один опрос
			pthread_mutex_lock(&mutex);
			bool isEmpty = tasks.empty();
			pthread_mutex_unlock(&mutex);
			//Если никто не дал задач
			if(isEmpty){
				#ifdef DEBUG
				printf("SENDING %d\n", m_rank);	
				#endif
				pthread_mutex_lock(&mutex);
				for(int p = 0; p < p_count; ++p){
					if(p != m_rank){
						MPI_Send(&PROC_OUT, 1, MPI_INT, p, TAG_GET, MPI_COMM_WORLD);
					}
				}
				workEnd = true;
				pthread_mutex_unlock(&mutex);

				#ifdef DEBUG
				printf("Observer end %d\n", m_rank);
				#endif
				break;
			}
		}
	}//while
}
	



/////////////////////////////////////////Observer


/////////////////////////////////////////Listener
void* startListen(void* t){
	pthread_mutex_lock(&mutex);
	if(p_count == 1 || lvl == 1){
		pthread_mutex_unlock(&mutex);

		#ifdef LVL_DEBUG
		printf("Listener lvl 1 end.\n");
		#endif

		return 0;
	} 
	pthread_mutex_unlock(&mutex);
	
	//deque<task>& tasks = *((deque<task>*)t);
	
	MPI_Status status;
	int count = 0;
	while(true){
		count++;
		int res = 0;
		MPI_Recv(&res, 1, MPI_INT, MPI_ANY_SOURCE, TAG_GET, MPI_COMM_WORLD, &status);


		//это будет верно только тогда, когда процесс
		//status.MPI_SOURCE отправит сообщение о завершении
		//по завершению всех поток закончит исполнение
		if(res == PROC_OUT){
			unAnsDied.erase(status.MPI_SOURCE);
			if(unAnsDied.empty()){
				#ifdef DEBUG
				printf("OUt: %d , %d\n",m_rank, count);
				#endif
				break;
			}
			//тк не надо отвечать,ведь прилетело из рассылки о завершении процесса
			//с рангом status.MPI_SOURCE
			continue;
		}


		pthread_mutex_lock(&mutex);
		//если завершил работу, но пришел запрос на задачи
		//отвечаем что задач нет и на другой стороне 
		//процесс уберет текущий из рассматриваемых на запрос
		if(workEnd){
			pthread_mutex_unlock(&mutex);
			MPI_Send(&PROC_OUT, 1, MPI_INT, status.MPI_SOURCE, TAG_ANS, MPI_COMM_WORLD);
			continue;
		}
		if(tasks.size() > res){
			task send_task = tasks.back();
			tasks.pop_back();
			pthread_mutex_unlock(&mutex);	
			MPI_Send(&send_task.time_mls, 1, MPI_INT, status.MPI_SOURCE, TAG_ANS, MPI_COMM_WORLD);
		}else{
			pthread_mutex_unlock(&mutex);	
			int count = 0;
			MPI_Send(&count, 1, MPI_INT, status.MPI_SOURCE, TAG_ANS, MPI_COMM_WORLD);
		}	
	}
	#ifdef DEBUG
	printf("Listener end %d\n", m_rank);
	#endif
}

/////////////////////////////////////////Listener



int main(int argc, char** argv)
{

	

	if(argc!=5){
		printf("Type tasks count than lvl than count of iterations.\n");
		return 0;
	}
	PROVIDED_TASKS = atoi(argv[1]);
	lvl = atoi(argv[2]);
	int type = atoi(argv[3]);
	int iteration =  atoi(argv[4]);

	int provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
	if(provided != MPI_THREAD_MULTIPLE){
		printf("WRONG INIT\n");
		return 0;
	}

	MPI_Comm_size(MPI_COMM_WORLD , &p_count);
	MPI_Comm_rank(MPI_COMM_WORLD , &m_rank);
	srand(m_rank);

	
	if(m_rank == 0)
		printf("Provided tasks:%d\n", PROVIDED_TASKS);
    
	//атрибуты потока
    pthread_attr_t attrs;
 
    //инициализация атрибутов потока
    if(0 != pthread_attr_init(&attrs))
    {
    	perror("Cannot initialize attributes");
    	abort();
    }

    if(0 != pthread_attr_setdetachstate(&attrs, PTHREAD_CREATE_JOINABLE))
    {
    	perror("Error in setting attributes");
   	 	abort();
    }
    
    
    std::chrono::time_point<std::chrono::system_clock> start, end;
    int glTime = 0;
	while(iteration--){
		for(int i = 0; i < p_count; ++i){
			proc.insert(i);
			unAnsDied.insert(i);
		}
		proc.erase(m_rank);
		unAnsDied.erase(m_rank);

		switch(type){
			case 0: initTasksLow(iteration);break;
			case 1: initTasksMedium(iteration);break;
			case 2: initTasksHigh(iteration);break;
		}

		//initTasksLow(iteration);
		//initTasksMedium(iteration);
		//initTasksHigh(iteration);


		//printTasks();

		start = std::chrono::system_clock::now();


		//MPI_Barrier(MPI_COMM_WORLD);
	    if(pthread_create(&worker, &attrs, startWork, 0)!=0)
	    {
	    	perror("Cannot create a thread worker");
	    	abort();
	    }

	    if(pthread_create(&listener, &attrs, startListen, 0) != 0)
	    {
	    	perror("Cannot create a thread listener");
	    	abort();
	    }

	    if(pthread_create(&observer, &attrs, startObserv, 0)!= 0)
	    {
	    	perror("Cannot create a thread listener");
	    	abort();
	    }

	    //pthread_attr_destroy(&attrs);

	    int w_ret = pthread_join(worker, NULL);
	    int l_ret = pthread_join(listener, NULL);
	    int o_ret = pthread_join(observer, NULL);
	    if(w_ret || l_ret || o_ret){
	    	perror("Cannot join a thread");
	    	abort();
	    }
	    //printf("WAIT\n");
	    //printf("rank:%d , tasks: %d\n",m_rank, finishedTasks);
	    MPI_Barrier(MPI_COMM_WORLD);


	    end = std::chrono::system_clock::now();
	    //seconds
	    int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>
                             (end-start).count();

	    //printf("P\n");
	    int done_tasks = 0;
	    int total_time = 0;

	    MPI_Reduce(&finishedTasks, &done_tasks, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	    MPI_Reduce(&elapsed, &total_time, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

	    if(m_rank == 0){
	    	if(done_tasks != PROVIDED_TASKS){
	    		//исправить завершение только 0
	    		printf("Miss tasks %d\n", done_tasks - PROVIDED_TASKS);
	    		break;
	    	}
	    	glTime+=elapsed;
	    	//printf("\nTotal Tasks: %d\n", done_tasks);
	    	//printf("Total Time: %d\n", total_time);
	    }


	    workEnd = false;
	    finishedTasks = 0;
	    globalWorkingTime = 0;
	    proc.clear();
	    if(!tasks.empty()){
	    	printf("TASK not empty\n");
	    	break;
	    }
	}
	if(m_rank == 0){
		printf("Global time: %d\n", glTime);
	}

	
    
    
 
    
    //printf("gone\n");
  
	MPI_Finalize();

    return 0;
}

