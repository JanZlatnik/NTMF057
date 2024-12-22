# Zápočtový program z NTMF057 Numerické metody pro teoretické fyziky I
## Úloha 2: vlastní vázané stavy radiálního problému

### Složka [`src`](src/)
  - modul [`core.rs`](src/core.rs) -> obsahuje základní fyzikální konstanty a inicializace nastavení
  - modul [`math.rs`](src/math.rs) -> zde jsou naprogramovány základní metody jako tvorba gridu, integrace, Runge-Kutte metoda, interpolace přirozeným kubickým splinem
  - modul [`radialSE.rs`](src/radialSE.rs) -> zde je metoda na řešení radiálního problému pomocí metody střelby
  - modul [`main.rs`](src/main.rs) -> špouští jednotlivé metody a ukládá data do souborů

### Aplikace
  - v hlavní složce je zkompolovaná aplikace [`Radial_problem.exe`](Radial_problem.exe), která ke spušnění vyžaduje přítomnost [`settings.toml`](settings.toml), ve kterém je uloženo nastavení, výstupy pak vypisuje do složek Task 1, Task 3 a Test
  - data ve složkách Task 1, Task 2 a Test byly spuštěny s příloženým nastavením

### složka [`Task 1`](Task%201/) obsahuje řešení první úlohy
  - [`RK_test.txt`](RK_test.txt) obsahuje test konvergence Runge-Kutteových metod
  - [`N2_potential.txt`](N2_potential.txt) obsahuje výpis potenciálu pro kontrolu
  - [`energies.txt`](energies.txt) obsahuje vlastní energie pro `l=0` a `l=10`
  - [`energies_control.txt`](energies_control.txt) obsahuje tytéž energie spočítané rozvojem do Fourierovy báze pro kontrolu
  - [`sign_function.txt`](sign_function.txt) obsahuje znaménko asymtotiky vlnové funkce propagované pomocí metody střelby jako kontrolní výstup při výpočtu pro `l=0`
  - [`eigenfunctions0.txt`](eigenfunctions0.txt) a [`eigenfunctions10.tx`](eigenfunctions10.txt) obsahují vlastní vlnové funkce pro `l=0` a `l=10`

· složka Task 3 obsahuje řešení třetí úlohy
  - interpolation.txt obsahuje potenciál a jeho aproximace přirozeným kubickým splinem s různým počtem bodů
  - interpolation_enegies.txt obsahuje vlastní energie pro l=0 spočítané pomocí metody střelby s přesným potenciálem a s potenciálem aproximovaným jednotlivými kubickými spliny
  - interpolation_energies_difference.txt obsahuje rozdíly vlastních energií pro l=0 spočítané pomocí metody střelby mezi přesným potenciálem a aproximací spliny
  - interpolation_ground_state_difference.txt obsahuje rozdíl energie základního stavu v závislosti na počtu použitých bodů pro interpolaci splinem

· složka Scripts obsahuje python skripty pro analýzu konvergence a tvorbu grafů
  - skript ConvergenceTest.py vytváří vždy složku s určitým nastavením a pak v ní spustí aplikaci
  - skripty ConvergenceGraphNR.py a ConvergenceGraphRmax.py z jednotlivých běhů sezbírají data z energies.txt pro l=0, porovnají je a vytvoří z nich grafy
  - skripty RungeKutta.py a spline.py pouze vytvářejí grafy ze souborů RK_test.txt a interpolation_ground_state_difference.txt

· složka Graphs obsahuje jednotlivé grafy
  - graf rmax_graph.pdf zobrazuje testování závislosti energií pro l=0 na koncovém bodu gridu, přičemž zorbazuje rozdíl mezi posledním během s rmax=12 au, počet bodů je vždy fixní na 1 au
    z grafu je vidět, že pro energie základního stavu se pro různé volby rmax liší až v řádu 10e-12, tedy s relativní chybou 10e-11
    pro vyšší stavy se rozdíly pro různé volby rmax pohybují mezi 10e-7 a 10e-9, tedy s relativní chybou 10e-5 až 10e-7
    pro několik nevyšších stavů můžeme vidět, že pro volby rmax 7.0 - 8.5 au nejsou výsledky ještě správně zkonvergované, neboť zde mají relativní chybu řádu 10e0, pro vyšší voplby rmax se již chovají podobně jako ostatní stavy
    => z hlediska volby rmax, pro rmax = 12 au, můžeme předpokládat, že energie vázaných stavů jsou určeny s relativní přesností zhrba 10e-7, u základního stavu dokonce 10e-11
  - graf nr_graph.pdf zobrazuje stejný test pro fixní rmax=12 au v závislosti na počtu bodů gridu
    => přesnost energií vázaných stavů je silně závislá na použitém počtu bodů grafu, přičemž můžeme říci, že řád metody by měl odpovídat řádu metody použité na řešení diferenciální rovnice - Runge-Kutta 4. řádu - neboť zvýšení počtu bodů o jeden řád odpovídá poklesu rozdílu energií o zhruba 4 řády
  - graf  RK_test.pdf zobrazuje testování Runge-Kutteových metod na řešení LHO na intervalu [0,20] s počátečními podmínkami nastavenými tak, aby přesné řešení odpovídalo A(t)=cos(t)
  - graf spline_interpolation.pdf zobrazuje vliv počtu podů použitých při interpolaci přirozenými kubickými spliny na rozdíl energií základního stavu
    => zde se nepřesnost také jeví jako O(1/N^4), což by odpovídalo nepřesnosti kubického splinu? Maximální počet použitých bodů pro spline je řádu 10e3, přičemž výpočet byl vždy prováděň na gridu s počtem bodů 10e4


