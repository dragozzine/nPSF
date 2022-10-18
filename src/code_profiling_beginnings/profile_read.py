# Cumulative knowledge document of commands used to profile code
# William Giforos, Summer 2022

import cProfile
import pstats
import io

p = pstats.Stats("profile.prof")

# .CUMULATIVE means we sort it according to cumulative time from longest to shortest
#p.print_stats();
#p.sort_stats(pstats.SortKey.CUMULATIVE).print_stats(3);
#p.sort_stats(pstats.SortKey.CUMULATIVE).print_stats(0.25);
#p.strip_dirs().sort_stats(pstats.SortKey.CUMULATIVE).print_stats(0.03);
#p.sort_stats('cumtime').print_stats(3);

#p.sort_stats('tottime').print_stats(3);

#p.strip_dirs().print_callers().sort_stats(pstats.SortKey.CUMULATIVE).print_stats(0.03);
#p.strip_dirs().sort_stats(pstats.SortKey.CUMULATIVE).print_stats(100);
p.strip_dirs().sort_stats(pstats.SortKey.TIME).print_stats(100);



#pr = cProfile.Profile()

##pr.dump_stats('profile1.prof')

#s1 = io.StringIO()
#ps1 = pstats.Stats(pr, stream=s1).sort_stats('cumtime')
#ps1.print_stats()

##with open('test_cum.txt', 'w+') as f:
 ##   f.write(s1.getvalue())

#s2 = io.StringIO()    
#ps2 = pstats.Stats(pr, stream=s2).sort_stats('tottime')
#ps2.print_stats()