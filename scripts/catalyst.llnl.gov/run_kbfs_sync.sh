#!/usr/bin/env bash

for n in 1 2 4 8 16 32
do
  rm -f sbatch${n}.out
  echo "#!/bin/bash" > run_${n}.sh
  echo "#SBATCH -N${n}" >> run_${n}.sh
  echo "#SBATCH -o sbatch${n}.out" >> run_${n}.sh
  echo "#SBATCH --ntasks-per-node=24" >> run_${n}.sh
  echo "#SBATCH -t 2:00:00" >> run_${n}.sh

  echo "export HAVOQGT_MAILBOX_SHM_SIZE=16384" >> run_${n}.sh
  echo "export HAVOQGT_MAILBOX_MPI_SIZE=131072" >> run_${n}.sh

  echo "srun --clear-ssd --ntasks-per-node=24 --distribution=block ./src/ingest_edge_list -o /dev/shm/graph -f 3 /p/lscratchf/havoqgtu/real/twitter/32x24parts/x0*" >> run_${n}.sh


  echo "export NUM_SOURCES=1" >> run_${n}.sh
  echo "srun -n1 -N1 sh scripts/do_cmake.sh" >> run_${n}.sh
  echo "srun -n1 -N1 make -j4" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 9178307" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 17481668" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 31181461" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 444682" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 29346643" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 19327398" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 16521058" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 26668632" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 11850563" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 20993397" >> run_${n}.sh


  echo "export NUM_SOURCES=2" >> run_${n}.sh
  echo "srun -n1 -N1 sh scripts/do_cmake.sh" >> run_${n}.sh
  echo "srun -n1 -N1 make -j4" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 6856398:20020227" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 16538908:5787088" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 16110355:10705414" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 25153063:12192027" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 5829602:352280" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 24736438:17554350" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 4631469:6672921" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 21475619:21236401" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 19759424:4341103" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 17279767:29446532" >> run_${n}.sh


  echo "export NUM_SOURCES=4" >> run_${n}.sh
  echo "srun -n1 -N1 sh scripts/do_cmake.sh" >> run_${n}.sh
  echo "srun -n1 -N1 make -j4" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 1138850:32425471:22427272:20588027" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 6358244:9927487:17836741:10088580" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 9249896:21228851:30519426:25991177" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 1679577:13850222:18411664:9392886" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 21185603:13882659:12578059:33125895" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 12155537:9533559:13949852:16187311" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 30194542:5317670:26243237:30327693" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 19743451:28902730:10753723:33440746" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 35224:28872311:11462463:28513387" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 15326060:15897284:33059149:22675817" >> run_${n}.sh


  echo "export NUM_SOURCES=8" >> run_${n}.sh
  echo "srun -n1 -N1 sh scripts/do_cmake.sh" >> run_${n}.sh
  echo "srun -n1 -N1 make -j4" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 8132757:32567175:19483970:33247779:12166052:30483589:17379962:4504829" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 10968283:1008119:4327654:4616503:20626840:1160403:31588445:331976" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 21029277:15033469:5823211:2003547:18060139:10105230:15163140:32005675" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 29138601:27082341:23772964:5823776:32766854:3201409:7492807:9889556" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 22622274:27008541:10072854:32558738:829526:29756384:23521924:6277410" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 4799382:20552129:25796355:30079396:8379473:23809283:970980:5922440" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 30455920:30669115:13867813:11365180:25216657:31435152:1701364:999060" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 6380931:29790334:14413520:20451152:3307733:4000314:532535:11637228" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 16213551:24396684:18751507:13439462:26712824:13555521:22497528:31575566" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 18625603:6245732:30133972:24631659:11373106:1074258:15016355:8789718" >> run_${n}.sh


  echo "export NUM_SOURCES=16" >> run_${n}.sh
  echo "srun -n1 -N1 sh scripts/do_cmake.sh" >> run_${n}.sh
  echo "srun -n1 -N1 make -j4" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 15733876:2180063:9580158:31255847:13226504:2018270:23956931:15034069:9257486:30957408:26441750:14456291:22171358:9169086:12953876:25943056" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 1728700:1654102:14422651:28012728:6238410:9634415:18533669:11778766:18110792:28436070:16915060:22759457:31350247:28172474:11306499:31715866" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 8603862:31682348:5303932:14174933:30504419:9003885:12076897:16734305:8403333:21198704:14123774:12760452:2992492:6667798:15022008:27977590" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 12335989:22700975:2339529:18053437:21056957:7542621:5007653:13151993:25851383:14649432:24079550:4298958:6802225:7434233:454592:18377527" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 8087018:6962887:4412847:11285275:11262811:14896871:833829:30192337:22930578:29753309:4175338:2173177:15313969:75148:5677850:22356742" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 721725:9499978:27822839:7892851:31948769:15584666:7411508:22362223:20263376:14791616:33506551:12292607:21690223:33151585:4422462:72310" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 599569:10304411:7507056:7697724:2449254:23195654:7380524:31069861:2424628:31239098:6250857:10089595:1642471:27058914:14266000:13106296" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 13295984:4005741:15757559:646857:33439858:23896349:17279113:13940253:25303826:14509381:30647141:32465863:24000076:13222392:15238976:25531102" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 6462157:13584026:20927070:28036167:23061776:25319703:17982189:17910679:19151168:3689371:24385808:9139206:7848006:5061633:5765424:31186085" >> run_${n}.sh
  echo "srun --drop-caches=pagecache -N${n} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s 22177545:24074315:14987744:25424492:5071836:20248241:11881676:20439268:24574548:31604672:30239410:17543906:6883741:3031960:11370565:18841207" >> run_${n}.sh

  sbatch run_${n}.sh

  sleep 5
done
