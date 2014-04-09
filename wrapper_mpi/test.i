#include"svipc.i"
shm_key=1234
sem_key=2345
write, "je suis dans Yorick!"
aaa=shm_read(shm_key,"donnee") 
pli, aaa;
for(i=1;i>0;i--){
  write, format="quit in %ds\n", i;
  pause, 1000;
}
aaa*=2;
shm_write,shm_key,"donnee", &aaa
pause, 100;
sem_give,sem_key,0;
quit;