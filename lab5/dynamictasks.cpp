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

using std::deque;
using std::set;
constexpr int PROC_OUT = -1;

//std::this_thread::sleep_for(std::chrono::milliseconds(x));

const int TAG = 12;
int p_count, m_rank;
pthread_mutex_t mutex;
pthread_t worker, listener, observer;
set<int> proc;


int globalWorkingTime = 0;//будет учитываться позже во время балансировки задач
int finishedTasks = 0;
bool workEnd = false , listenerEnd = false;


typedef struct Task{
	int time_mls;
}task;



/////////////////////////////////////////INIT
void initTasks(deque<task>& tasks){
	task t;
	for (int i = 0; i < m_rank + 5; ++i)
	{
		t.time_mls = 1 + 100 * m_rank * 2;
		tasks.push_back(t);
	}
}

void printTasks(deque<task>& tasks){
	for (int i = 0; i < tasks.size(); ++i)
	{
		//printf("Task %d: %d mls \n", i, tasks[i].time_mls);
	}
}
/////////////////////////////////////////INIT



/////////////////////////////////////////Worker


void* startWork(void* t){
	deque<task>& tasks = *((deque<task>*)t);

	while(true){

		pthread_mutex_lock(&mutex);
		if(workEnd){
			pthread_mutex_unlock(&mutex);
			break;
		}
		pthread_mutex_unlock(&mutex);

		task* cur_task = nullptr;
		pthread_mutex_lock(&mutex);
		//printf("%d: %d\n", m_rank, (int)tasks.size());
		if(!tasks.empty()){
			cur_task = &tasks.front();
			tasks.pop_front();
		}
		pthread_mutex_unlock(&mutex);
		if(cur_task != nullptr){
			std::this_thread::sleep_for(std::chrono::milliseconds(cur_task->time_mls));

			pthread_mutex_lock(&mutex);
			globalWorkingTime+=cur_task->time_mls;
			pthread_mutex_unlock(&mutex);
			++finishedTasks;
		}
	}
	//printf("%d , Worker end\n", m_rank);
}




/////////////////////////////////////////Worker



/////////////////////////////////////////Observer




void* startObserv(void* t){
	deque<task>& tasks = *((deque<task>*)t);
	MPI_Status status;
	set<int>::iterator it;
	while(true){
		int tasks_count = 0;
		pthread_mutex_lock(&mutex);
		if(workEnd){
			pthread_mutex_unlock(&mutex);
			break;
		}
		tasks_count = tasks.size();
		pthread_mutex_unlock(&mutex);

		if(tasks_count < 2){
			for (it=proc.begin(); it!=proc.end(); ++it){

    			int i = *it;
				int res = 0;
				//будет передавать кол-во ост задач



				MPI_Send(&tasks_count, 1, MPI_INT, i, 12, MPI_COMM_WORLD);
				MPI_Recv(&res, 1, MPI_INT, i, 21, MPI_COMM_WORLD, &status);
				
				if(res > 0){
					task new_task;
					new_task.time_mls = res;

					pthread_mutex_lock(&mutex);
					globalWorkingTime+=res;
					tasks.push_back(new_task);
					pthread_mutex_unlock(&mutex);
					break;
				}else if(res == PROC_OUT){
					//printf("Recv OUT\n");
					proc.erase(i);
					//it++;
				}

			}//for
				//за время опроса могли кончиться задачи
				//поэтому проверим что их нет 
				//для того чтобы не проводить еще один опрос
				pthread_mutex_lock(&mutex);
				tasks_count = tasks.size();
				pthread_mutex_unlock(&mutex);
				//Если никто не дал задач
				if(tasks_count == 0){					
					for (int p =0; p < p_count; ++p){
						if(p != m_rank)
							MPI_Send(&PROC_OUT, 1, MPI_INT, p, 12, MPI_COMM_WORLD);
					}
					pthread_mutex_lock(&mutex);
					workEnd = true;
					pthread_mutex_unlock(&mutex);

					//printf("Observer end %d\n", m_rank);
					break;
				}
		}
	}//while
}
	



/////////////////////////////////////////Observer


/////////////////////////////////////////Listener
void* startListen(void* t){
	deque<task>& tasks = *((deque<task>*)t);
	
	MPI_Status status;

	while(true){

		int res = 0;
		MPI_Recv(&res, 1, MPI_INT, MPI_ANY_SOURCE, 12, MPI_COMM_WORLD, &status);

		if(res == PROC_OUT){
			pthread_mutex_lock(&mutex);
			proc.erase(status.MPI_SOURCE);
			bool b = proc.empty();
			pthread_mutex_unlock(&mutex);
			if(b)break;
		}

		pthread_mutex_lock(&mutex);
		if(workEnd){
			//printf("OUT\n");
			pthread_mutex_unlock(&mutex);
			MPI_Send(&PROC_OUT, 1, MPI_INT, status.MPI_SOURCE, 21, MPI_COMM_WORLD);
			proc.erase(status.MPI_SOURCE);
			if(proc.empty()){
				pthread_mutex_unlock(&mutex);
				 break;
			}
		}

		
		if(tasks.size() > res){
			//можно сделать проверку на то, чтобы не посылать задачку, если она одна
			task send_task = tasks.back();
			tasks.pop_back();
			pthread_mutex_unlock(&mutex);	
			MPI_Send(&send_task.time_mls, 1, MPI_INT, status.MPI_SOURCE, 21, MPI_COMM_WORLD);
		}else{
			pthread_mutex_unlock(&mutex);	
			int count = 0;
			MPI_Send(&count, 1, MPI_INT, status.MPI_SOURCE, 21, MPI_COMM_WORLD);
		}	
	}
	//printf("Listener end %d\n", m_rank);
}

/////////////////////////////////////////Listener



int main(int argc, char** argv)
{
	//MPI_Init(&argc  , &argv);
	int provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
	if(provided != MPI_THREAD_MULTIPLE)return 0;

	MPI_Comm_size(MPI_COMM_WORLD , &p_count);
	MPI_Comm_rank(MPI_COMM_WORLD , &m_rank);

	for(int i = 0; i < p_count; ++i){
		proc.insert(i);
	}
	proc.erase(m_rank);
	//set<int>::iterator it;

	/*for(it=proc.begin(); it!=proc.end(); ++it){
		printf("%d\n", *it);
		if(*it == 2){
			proc.erase(0);
			proc.erase(1);
			proc.erase(2);
		}
	}*/


	deque<task> tasks;

	initTasks(tasks);
	printTasks(tasks);


	//std::this_thread::sleep_for(std::chrono::milliseconds(1000));

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

    if(pthread_create(&worker, &attrs, startWork, (void*)&tasks) != 0)
    {
    	perror("Cannot create a thread worker");
    	abort();
    }

    if(pthread_create(&listener, &attrs, startListen, (void*)&tasks) != 0)
    {
    	perror("Cannot create a thread listener");
    	abort();
    }

    if(pthread_create(&observer, &attrs, startObserv, (void*)&tasks) != 0)
    {
    	perror("Cannot create a thread listener");
    	abort();
    }

    pthread_attr_destroy(&attrs);

    int w_ret = pthread_join(worker, NULL);
    int l_ret = pthread_join(listener, NULL);
    int o_ret = pthread_join(observer, NULL);
    if(w_ret || l_ret || o_ret){
    	perror("Cannot join a thread");
    	abort();
    }

    printf("\n%d: %d\n", m_rank, finishedTasks);


	MPI_Finalize();
    return 0;
}

