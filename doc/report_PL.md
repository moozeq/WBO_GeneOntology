# Raport

Repozytorium wraz z poniższym raportem dostępne jest również pod [tym linkiem](https://github.com/moozeq/WBO_GeneOntology).

## Część A

### 1

Aby wykonać wyszukiwanie najbliższych sekwencji do zadanych fragmentów białek, użyta została
lokalna wersja programu `BLAST`, która dla sekwencji nukleotydowych genów _E. coli_ zbudowała bazę oraz
skrypt [hammer.py](hammer.py), który wywołuje program `BLAST`.
Do tworzenia bazy i wyszukiwania sekwencji użyty został interfejs `Bio.Blast.Applications` z biblioteki
`Biopython`. Do wyszukiwania sekwencji użyto wersji `tblastn` ponieważ fragmenty były sekwencjami
aminokwasowymi, a baza stworzona została dla sekwencji nukleotydowych. Głównym parametrem przy
przeszukiwaniu było `e-value`, które ustalone zostało na `1e-20` ponieważ zależało nam
na znalezieniu najbardziej prawdopodobnych sekwencji (współczynnik ustalony został na podstawie
[odpowiedzi na podobne pytanie na stackexchange](https://biology.stackexchange.com/questions/19338/what-e-value-blast-cut-off-is-best),
oraz manualnym przeglądaniu wyników bez obcięcia). Ponieważ przy zbieraniu wyników zależało nam
na otrzymaniu genów, do których pasowały zadane fragmenty, został również określony parametr
identyczności na 90%, z tych samych powodów co powyżej. Wyniki analizy zapisywane są do pliku `results/blasted.fasta`.

### 2

Wyszukiwanie domen PFAM zostało przeprowadzone przy użyciu lokalnej wersji programu `HMMER` i skryptu
[blaster.py](blaster.py).
Do wywołań programu `HMMER` użyte zostały funkcje z biblioteki `subprocess`, do parsowania
pliku wynikowego użyty został interfejs `Bio.SearchIO.HmmerIO` z biblioteki `Biopython`.
Tak samo jak przy użyciu aplikacji `BLAST` - `e-value` został określony na `1e-20` przy ręcznym
przeglądaniu wyników. Dzięki temu otrzymano listę genów, w których odnaleziono domeny PFAM wraz z adnotacjami terminów.
Wyniki analizy zapisywane są do pliku `results/study_annotated.json` (znalezione geny wraz z adnotacjami
z programu `HMMER`) oraz `results/studyPfam.txt` (znalezione geny zawierające domeny PFAM do dalszej analizy).

### 3

Poszukiwanie terminów zostało osiągnięte dzięki skryptowi [godb.py](godb.py) oraz bibliotece [goatools](https://github.com/tanghaibao/goatools).
Skrypt odpowiada za pobranie plik bazy `go.obo` do wyszukania terminów na postawie identyfikatora GO, a następnie
na podstawie otrzymanych danych z biblioteki `goatools` (mapy genów na zestawy terminów GO `results/ecoli.assocs`), odpowiedniej adnotacji znalezionych genów co
zakończone jest modyfikacją pliku `results/study_annotated.json` zawierającego znalezione geny z adnotacjami
zarówno z programu `HMMER` jak i z bazy `GO`.

### 4

Analiza nadreprezentacji genów ze zbiorów `A`, `B` oraz znalezionych w powyższej analizie genów zostało
wykonane przy użyciu biblioteki `goatools`, która jest bardzo prosta w użyciu, a daje nam wyniki używając
testu Fisher'a z poprawką Bonferroniego, co było wymaganiem w projekcie. Parametrem, który należało określić
było tu `p-value`, które pozostało na domyślnej wartości proponowanej przez bibliotekę `goatools`: `0.05`.
Do podzielenia genów na odpowiednie grupy został użyty skrypt [splitter.py](splitter.py), dzięki któremu
otrzymane zostały 3 pliki:
- `results/studyA.txt`: lista genów z grupy `A`
- `results/studyB.txt`: lista genów z grupy `B`
- `results/population.txt` lista wszystkich genów

Dołączony został do nich również plik `results/studyPfam.txt` zawierający nazwy genów znalezionych w poprzednich 
krokach analizy.

Wyniki analizy nadreprezentacji terminów adnotacji są zapisywane do następujących plików:
- `results/enrichmentA.tsv`: analiza nadreprezentacji dla genów z grupy `A`
- `results/enrichmentB.tsv`: analiza nadreprezentacji dla genów z grupy `B`
- `results/enrichmentPfam.tsv` analiza nadreprezentacji dla genów znalezionych w poprzednich krokach analizy

## Część B

### 1

Przy wyszukiwaniu sekwencji otrzymano plik `results/blasted.tsv`, który zawierał wszystkie znalezione trafienia
w ilości `157`. Po przefiltrowaniu go względem identyczności znalezione geny zostały zapisane wraz
z ich sekwencjami do pliku `results/blasted.fasta`, który zawierał równo `100` sekwencji aminokwasowych.

### 2

Przy przeszukiwaniu `100` sekwencji genów pod kątem zawierania prawdopodobnych domen PFAM z `e-value`
określonym na `1e-20` znaleziono `132` domeny PFAM w `90` z nich, wyniki wraz z adnotacjami zostały zapisane do
pliku `results/FINAL_STUDY_ANNOTATED.json`.

### 3

Korzystając z pliku `results/ecoli.assocs` do `90` genów zawierających prawdopodobne domeny PFAM, przypisano równo
`1000` rozwiniętych terminów `GO` wraz z ich identyfikatorami.

### 4

Dla białek z grup `A`, `B` i znalezionych na etapach wcześniejszej analizy wykonane zostały analizy nadreprezentacji
terminów ich funkcji używając testu Fisher'a (licząc przy tym także poprawki m.in. Bonferroniego). Przy `p-value` równym `0.05` znaleziono
następującą liczbę nadreprezentacji terminów (zebrane wyniki znajdują się w pliku `results/FINAL_ENRICHMENT.json`):

|grupa badana|grupa kontrolna|liczba nadreprezentowanych terminów|
|:---:|:---:|:---:|
|A|populacja|145|
|B|populacja|91|
|Pfam|populacja|186|
|A|B|21|

Po posortowaniu terminów względem najniższych `p-value` i wzięciu pierwszych do 5 wyników, dla każdej z kategorii
(BP=Biological Process; MF=Molecular Function; CC=Cellular Component), otrzymujemy poniższe wyniki (dostępne także w plikach `results/FINAL_ENRICHMENT_*.tsv`):

_Z poprawką Bonferroniego dla `p-value`
równego `0.5` udało się znaleźć jedynie bardzo ogólne właściwości, np. że produkty genów są komponentem komórki, czy, że znajdują się w cytozolu._

#### Wzbogacenie w genach grupy A do populacji

|GO|kategoria|opis|liczba w badaniu|liczba w próbie kontrolnej|p-value|
|:---:|:---:|:---:|:---:|:---:|:---:|
|GO:0097428|BP|protein maturation by iron-sulfur cluster transfer|5/466|6/2550|0.0010201587697092848|
|GO:0046471|BP|phosphatidylglycerol metabolic process|4/466|5/2550|0.00471682516496374|
|GO:0032465|BP|regulation of cytokinesis|3/466|3/2550|0.006070820759740798|
|GO:0046486|BP|glycerolipid metabolic process|7/466|14/2550|0.006832865597283435|
|GO:0006650|BP|glycerophospholipid metabolic process|7/466|14/2550|0.006832865597283435|
|GO:0051287|MF|NAD binding|17/466|44/2550|0.0012399890059333499|
|GO:0003824|MF|catalytic activity|282/466|1378/2550|0.002024038358785938|
|GO:0070403|MF|NAD+ binding|5/466|7/2550|0.003034144470912581|
|GO:0042802|MF|identical protein binding|84/466|349/2550|0.003561591603824005|
|GO:0035251|MF|UDP-glucosyltransferase activity|3/466|3/2550|0.006070820759740798|
|GO:0110165|CC|cellular anatomical entity|375/466|1887/2550|0.00035954154443360785|
|GO:0005575|CC|cellular_component|384/466|1948/2550|0.0005801999624080763|
|GO:0005829|CC|cytosol|168/466|770/2550|0.0025830100279591885|
|GO:0005737|CC|cytoplasm|155/466|702/2550|0.002817150640428286|
|GO:0043228|CC|non-membrane-bounded organelle|7/466|96/2550|0.0028417708105327793|

Z powyższej analizy możemy przypuszczać, iż produkty genów grupy `A` znajdują się w cytozolu komórki,
pełnią funkcję katalityczną, związaną z procesami metabolitycznymi odpowiedzialnymi za zarządzanie
energią w komórce.

Z poprawką Bonferroniego:

|GO|kategoria|opis|liczba w badaniu|liczba w próbie kontrolnej|p-value|p-value Bonferroni|
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
|GO:0110165|CC|cellular anatomical entity|375/466|1887/2550|0.00035954154443360785|0.0913235522861364|
|GO:0005575|CC|cellular_component|384/466|1948/2550|0.0005801999624080763|0.14737079045165138|

#### Wzbogacenie w genach grupy B do populacji

|GO|kategoria|opis|liczba w badaniu|liczba w próbie kontrolnej|p-value|
|:---:|:---:|:---:|:---:|:---:|:---:|
|GO:2000113|BP|negative regulation of cellular macromolecule biosynthetic process|3/336|107/2550|0.00035504187807809337|
|GO:0051171|BP|regulation of nitrogen compound metabolic process|20/336|288/2550|0.0005730487078432635|
|GO:0010558|BP|negative regulation of macromolecule biosynthetic process|4/336|114/2550|0.0005938337036435375|
|GO:0031327|BP|negative regulation of cellular biosynthetic process|4/336|115/2550|0.0005940861422437907|
|GO:0009890|BP|negative regulation of biosynthetic process|4/336|115/2550|0.0005940861422437907|
|GO:0140110|MF|transcription regulator activity|7/336|165/2550|0.00012080546433101622|
|GO:0003700|MF|DNA-binding transcription factor activity|7/336|147/2550|0.0009562953601487175|
|GO:0003824|MF|catalytic activity|209/336|1378/2550|0.0014857153789797491|
|GO:0005353|MF|fructose transmembrane transporter activity|4/336|7/2550|0.007475779631354398|
|GO:0022877|MF|protein-N(PI)-phosphohistidine-fructose phosphotransferase system transporter activity|4/336|7/2550|0.007475779631354398|
|GO:0005575|CC|cellular_component|284/336|1948/2550|0.00010845204565379978|
|GO:0110165|CC|cellular anatomical entity|276/336|1887/2550|0.0002300552060598511|
|GO:0005829|CC|cytosol|122/336|770/2550|0.010656334636478446|
|GO:0031225|CC|anchored component of membrane|4/336|10/2550|0.032444079876899554

Z powyższej analizy możemy przypuszczać, iż produkty genów grupy `B` znajdują się także w cytozolu komórki, bądź
zakotwiczone są w membranie, wiążą się do DNA i związane są z procesami regulacji syntezy produktów oraz
regulacji transkrypcji.

Z poprawką Bonferroniego:

|GO|kategoria|opis|liczba w badaniu|liczba w próbie kontrolnej|p-value|p-value Bonferroni|
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
|GO:0140110|MF|transcription regulator activity|7/336|165/2550|0.00012080546433101622|0.26589282699256667|
|GO:0005575|CC|cellular_component|284/336|1948/2550|0.00010845204565379978|0.027546819596065142|
|GO:0110165|CC|cellular anatomical entity|276/336|1887/2550|0.0002300552060598511|0.05843402233920218|


#### Porównanie wzbogacenia pomiędzy zbiorami A i B

|GO|kategoria|opis|liczba w badaniu|liczba w próbie kontrolnej|p-value|
|:---:|:---:|:---:|:---:|:---:|:---:|
|GO:0006979|BP|response to oxidative stress|15/466|17/802|0.011499487157127465|
|GO:0010556|BP|regulation of macromolecule biosynthetic process|47/466|65/802|0.017826485389559502|
|GO:2000113|BP|negative regulation of cellular macromolecule biosynthetic process|16/466|19/802|0.019187132645730666|
|GO:2000112|BP|regulation of cellular macromolecule biosynthetic process|44/466|61/802|0.021787216273690478|
|GO:0051171|BP|regulation of nitrogen compound metabolic process|50/466|70/802|0.021948726462711978|
|GO:0000287|MF|magnesium ion binding|19/466|44/802|0.042052687751242934

Z powyższej analizy możemy przypuszczać, iż wspólną funkcją pełnioną przez produkty genów z grup `A` i `B`
jest regulacja biosyntezy niektórych składników komórki i odpowiedzi na stres oksydacyjny, przy czym
pełnią funkcję składnika wiążącego jony magnezu.

Z poprawką Bonferroniego nie udało się uzyskać w tym wypadku żadnych wyników.

#### Wzbogacenie w genach znalezionych w analizie do populacji

|GO|kategoria|opis|liczba w badaniu|liczba w próbie kontrolnej|p-value|
|:---:|:---:|:---:|:---:|:---:|:---:|
|GO:0071103|BP|DNA conformation change|6/90|25/2550|0.00016857988976947191|
|GO:0051276|BP|chromosome organization|6/90|30/2550|0.0004900889563211258|
|GO:0006323|BP|DNA packaging|3/90|6/2550|0.0007875398115891808|
|GO:0006265|BP|DNA topological change|3/90|6/2550|0.0007875398115891808|
|GO:0030261|BP|chromosome condensation|3/90|6/2550|0.0007875398115891808|
|GO:0016491|MF|oxidoreductase activity|26/90|318/2550|2.2892324880787546e-05|
|GO:0003674|MF|molecular_function|89/90|2198/2550|3.2213492365591133e-05|
|GO:0097159|MF|organic cyclic compound binding|51/90|993/2550|0.0005916419422256992|
|GO:1901363|MF|heterocyclic compound binding|51/90|993/2550|0.0005916419422256992|
|GO:0003916|MF|DNA topoisomerase activity|3/90|6/2550|0.0007875398115891808|
|GO:0005829|CC|cytosol|42/90|770/2550|0.0009474822386315042|
|GO:0005694|CC|chromosome|3/90|7/2550|0.0013430415733970092|
|GO:1990204|CC|oxidoreductase complex|8/90|63/2550|0.0013854772560385157|
|GO:0009330|CC|DNA topoisomerase type II (double strand cut, ATP-hydrolyzing) complex|2/90|3/2550|0.003611830705684397|
|GO:0009295|CC|nucleoid|3/90|17/2550|0.02019631865851119|

Z powyższej analizy możemy przypuszczać, iż badane przez nas fragmenty białek, które były produktami
odpowiednich genów w komórce, odpowiadają głównie za regulację konformacji DNA i organizację chromosomów w komórce.

Z poprawką Bonferroniego:

|GO|kategoria|opis|liczba w badaniu|liczba w próbie kontrolnej|p-value|p-value Bonferroni|
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
|GO:0071103|BP|DNA conformation change|6/90|25/2550|0.00016857988976947191|0.48770162110308224|
|GO:0016491|MF|oxidoreductase activity|26/90|318/2550|2.2892324880787546e-05|0.05038600706261339|
|GO:0003674|MF|molecular_function|89/90|2198/2550|3.2213492365591133e-05|0.07090189669666608|
|GO:0005829|CC|cytosol|42/90|770/2550|0.0009474822386315042|0.24066048861240208|
|GO:0005694|CC|chromosome|3/90|7/2550|0.0013430415733970092|0.34113255964284034|
|GO:1990204|CC|oxidoreductase complex|8/90|63/2550|0.0013854772560385157|0.35191122303378297|
