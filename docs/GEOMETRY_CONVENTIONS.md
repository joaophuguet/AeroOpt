# Convenções geométricas — v0.1.1

## Referenciais e unidades

O referencial global tem x positivo a jusante e y positivo para cima. A corda
principal é a unidade: `c_main = 1`. A referência de ataque principal é `(0,0)`
e o ponto médio de fuga principal é `(1,0)`. A palavra referência é intencional:
o contorno amostrado não define a curva contínua exata do nariz.

Coordenadas, posição nominal, deslocamentos e distância mínima são expressos
em `c_main`; áreas globais, em `c_main²`. Ângulos da API e do JSON são graus.
`chord_ratio = c_f/c_main` é adimensional e estritamente positivo. Não há
conversão implícita para metros, milímetros ou radianos. Cada arquivo original
deve usar uma única unidade de comprimento, igual nos dois eixos. Os dois
arquivos podem ter escalas e números de pontos diferentes: cada qual recebe
sua própria normalização antes da montagem.

## Contrato do contorno e referências

Um contorno é uma sequência de pontos em ordem de superfície, começando em um
lado da borda de fuga, passando pelo nariz e terminando no outro lado da fuga.
Pode-se percorrer primeiro o lado superior ou o inferior. A primeira e a última
linha numérica definem as duas extremidades de fuga. Seu ponto médio é

`T = p_first + (p_last - p_first)/2`.

Uma fuga de espessura nula pode repetir exatamente o primeiro ponto na última
linha. Essa única duplicação é permitida e permanece nos arrays. Nenhuma outra
duplicação ou quase duplicação é removida automaticamente. Uma fuga de espessura
finita possui duas extremidades diferentes.

Por padrão, a referência de ataque é o **ponto amostrado mais distante de T**,
usando distância euclidiana, e não a menor coordenada x. Em empate, vence o
primeiro índice (comportamento determinístico de `argmax`). O índice precisa
estar no interior da sequência. Esse procedimento é invariante a translação,
rotação própria e escala uniforme, salvo arredondamento, e lida com pequenas
imperfeições sem impor uma orientação inicial do arquivo.

Esse critério **não identifica uma LE geométrica exata para qualquer forma**.
Uma distribuição assimétrica, nariz pouco amostrado, camber extremo ou forma
incomum pode deslocar o ponto mais distante para fora da referência desejada.
Não se ajusta círculo, spline nem linha média. O usuário pode escolher
`Airfoil.from_dat(..., leading_edge_index=i)` ou os argumentos
`main_leading_edge_index` e `flap_leading_edge_index` em `from_dat` da montagem.
O índice é base zero nas linhas numéricas, ignorando cabeçalho e linhas vazias.
O método e o índice efetivos são registrados no JSON.

O programa não pode deduzir semanticamente a fuga de um contorno arbitrário:
uma sequência que começa na LE ou no meio de uma superfície viola o contrato,
mesmo que gere um polígono válido. Não há reordenação automática ou garantia
de detecção de toda violação. Antes de usar novos conjuntos de perfis, revise
as extremidades e referências. Um par inicial de inteiros positivos ambos >1
é rejeitado no parser por ambiguidade com contagens Lednicer; use o construtor
explícito por pontos quando essa coincidência for uma coordenada legítima, ou
declare `AeroOpt_COORDINATES_v0.1` como primeiro token do cabeçalho. Esse marcador
é uma declaração explícita de formato, não uma tentativa de inferir contagens.
O exportador sempre o emite para permitir round trip de posições globais
arbitrárias, inclusive quando a primeira coordenada é um par de inteiros >1.

## Normalização

Seja `L = p[i_LE]`, `c = ||T-L||`, `e_x = (T-L)/c` e
`e_y = (-e_x_y, e_x_x)`. Para cada ponto original:

`q = [ dot(p-L, e_x), dot(p-L, e_y) ] / c`.

Essa transformação rígida própria seguida de escala uniforme leva L a `(0,0)`
e T a `(1,0)`. Distâncias ficam divididas por c e áreas por c²; ângulos e forma
são preservados. A transformação não espelha o perfil. Inverter a ordem de
percurso não muda o lado físico de cada ponto. Reflexões não são corrigidas.

Os pontos individuais de uma fuga aberta não precisam ter x exatamente 1:
quem ocupa `(1,0)` é seu ponto médio. Da mesma forma, uma imperfeição pode deixar
outro ponto com x ligeiramente negativo. Não se recortam coordenadas nem se
forçam pontos extremos. Não há suavização, interpolação ou redistribuição.

`original_points` guarda os números lidos e `normalized_points` um array
separado. Ambos são copiados para buffers somente leitura e protegidos contra
reativação de escrita. A corda, L e T permanecem disponíveis no referencial
original. Nomes e caminhos absolutos são preservados como proveniência.

## Colocação do flap

Após normalizar o flap na própria corda, `q_LE=(0,0)`. Seja
`r = c_f/c_main`, `a=(x_nom,y_nom)`, `d=(dx,dy)` e `delta` em graus.

1. Escala uniforme: `s = r*q` (agora em cordas principais).
2. Rotação em torno da LE local: `s_rot = R(-delta)*s`.
3. Translação: `p_global = a + d + s_rot`.

Para vetores coluna:

```text
R(-delta) = [ cos(delta)   sin(delta) ]
            [-sin(delta)   cos(delta) ]
```

A implementação usa arrays de vetores linha e multiplica por `R.T`. Deflexão
positiva é horária: o ponto local de fuga `(1,0)`, com `delta=90°`, fica abaixo
da LE em `(0,-r)` antes da translação. A posição efetiva da LE é sempre `a+d`.
Os deslocamentos são aplicados nos eixos globais, não nos eixos girados do flap.
Essa é a fórmula solicitada com `p_ref_local=(0,0)` na geometria escalonada e
`p_ref_global=a+d`. O ângulo é reduzido módulo 360 apenas para a avaliação
trigonométrica; o parâmetro original é preservado nos metadados.

`configure` cria outra configuração com valores absolutos, recalculada a partir
de `flap.normalized_points`. Não gira uma configuração anterior. Argumentos
omitidos voltam a zero. A base não é alterada. Todos os parâmetros são finitos;
escala e transformações não representáveis são rejeitadas.

## Polígonos sólidos, fuga aberta e áreas

Para Shapely, área, interseção e distância, os pontos definem uma poligonal
fechada por **um segmento reto entre o último e o primeiro ponto**. Uma fuga
aberta recebe, assim, uma tampa geométrica de espessura finita. Essa tampa
participa da área, da distância e das interseções. Não representa uma hipótese
de escoamento. Os arquivos exportados mantêm as linhas originais transformadas,
sem adicionar amostras; a mesma regra de fechamento é declarada nos metadados
e aplicada nas figuras.

Área é a área não orientada do polígono válido, nunca uma soma algébrica de
lóbulos autointersectantes. `Airfoil.area` e `normalized_area` retornam `None`
quando o polígono é inválido ou tem área desprezível. Um contorno pode ter
envoltória não degenerada, mas autointersectar: ele é mantido para diagnóstico,
recebe rejeição estruturada em `validate` e não pode ser exportado como válido.

## Validação, distância e tolerâncias

Importação/construção rejeita menos de três pares, valores não finitos,
duplicações problemáticas, envoltória colinear ou desprezível e corda degenerada.
As verificações iniciais de proximidade usam distâncias relativas à diagonal
da caixa dos dados: `LENGTH_TOL=1e-10`; a área da envoltória nesse referencial
deve exceder `AREA_TOL=1e-12`. Isso evita que uma simples mudança de escala
original determine a aceitação. A exceção para duplicação final exige igualdade
exata, sem substituir silenciosamente um ponto quase coincidente.

Na v0.1.1 a montagem distingue duas escalas. Qualidade intrínseca de cada
contorno usa `D = diagonal da própria caixa delimitadora`: duplicações são
verificadas em `(p-p_first)/D`, com `LENGTH_TOL=1e-10`, e a área nesse referencial
deve exceder `AREA_TOL=1e-12`. Ou seja, os limites equivalem a `1e-10 D` e
`1e-12 D²`. A área de validação é calculada nesse referencial local para não
subtrair números globais grandes desnecessariamente. Polígonos usados nas
operações entre elementos continuam com os pontos globais originais.

Contato e folga usam **sempre** `LENGTH_TOL=1e-10 c_main`. As constantes não
mudaram; corrigiu-se sua aplicação à qualidade intrínseca do elemento. Na v0.1,
uma escala r multiplicava a área por r² mas deixava fixa a área mínima global,
gerando rejeição artificial; pontos legítimos também eram classificados como
duplicados quando o elemento era reduzido. Um contorno resolvido deve manter
seu diagnóstico intrínseco sob escala uniforme, pois `A/D²` e `distância/D` são
invariantes nessa operação. Há testes reproduzindo ambos os defeitos
antes/depois da correção. A diagonal da caixa é uma escala numérica conveniente,
não um diâmetro invariante a rotação; casos intrínsecos exatamente próximos
dos limiares exigem revisão, especialmente em contornos muito pouco resolvidos.

Essas políticas não são incertezas físicas. Coordenadas globais que perdem
resolução em float64 ainda podem colapsar pontos e ser rejeitadas corretamente:
nenhum reescalonamento recupera informação já perdida. O teste com razão de
corda `1e-7` é uma sondagem numérica, não uma recomendação de projeto. Contornos
repetidos, área relativa muito pequena e autointerseção são diagnosticados
antes de operar sobre pares de sólidos inválidos. `ValidationReport` mantém
os campos existentes, acrescentando unidades explícitas para cada tolerância.

A distância mínima global é `distance(boundary(main), boundary(flap))`, incluindo
as tampas de fuga. É não negativa e está definida apenas quando os dois
elementos são válidos isoladamente. Para sólidos aninhados, a distância entre
contornos pode ser positiva, **apesar de a configuração ser inválida por
sobreposição sólida**. Por isso a distância não substitui a operação de
interseção dos polígonos.

Contato ou cruzamento exato de contornos é rejeitado mesmo com folga solicitada
zero. Toda área de interseção estritamente positiva é rejeitada, inclusive
contenção; não se usa a tolerância de área para autorizar pequenas penetrações.
Separação disjunta <= `LENGTH_TOL` também é rejeitada como `near_contact`.
Uma folga configurável `h >= 0` é insuficiente se `distance + LENGTH_TOL < h`.
Assim, igualdade é aceita dentro dessa tolerância; isso precisa ser considerado
ao escolher margens físicas. Não há `buffer(0)`, snap, reparo ou arredondamento
de contornos. Consulte `ValidationReport.issues`, não apenas `valid`.

## Gap e overlap aerodinâmicos

A v0.1 não define gap nem overlap aerodinâmicos. Referências de fuga principal
e ataque do flap estão disponíveis, mas sua separação em x/y não determina
por si só uma garganta de escoamento ou uma convenção aerodinâmica universal.
Rotação, espessura de fuga e formas arbitrárias podem mudar a região relevante.

Esses campos são `null` no JSON, com `gap_overlap_status=pending_definition`.
`overlap_area` é apenas a área de interseção **sólida**, em cordas principais ao
quadrado; não é overlap aerodinâmico. Antes de adicionar tais métricas, definir
e revisar referências, direção de medição, região de busca, sinais, unidades,
casos de falha e compatibilidade com a literatura/solver escolhido.

## Exportação e reprodutibilidade

DAT separado por elemento, ambos no referencial global, com 17 algarismos
significativos; JSON UTF-8 sem NaN com esquema `aeroopt.geometry.v0.1`.
Parâmetros, identificação, origem, contagens, método de normalização,
referências, unidades, versões de Python/GEOS/bibliotecas e relatório acompanham a montagem.
O array serializado é o mesmo usado na validação da exportação. A ordem dos
pontos não muda. Reimportar o DAT em `original_points` recupera as coordenadas
globais; pedir `normalized_points` executa intencionalmente outra normalização.

Não há números aleatórios nem timestamps nos DAT/JSON de uma configuração.
Mesmos perfis, parâmetros, identificação e versões produzem os mesmos dados;
o caminho absoluto de origem pode mudar os metadados entre máquinas. Não se
promete identidade byte a byte de figuras SVG/PNG (metadados do backend/fontes)
nem igualdade bit a bit entre diferentes versões de GEOS, NumPy ou plataformas.

## Referências de implementação

- [Shapely User Manual](https://shapely.readthedocs.io/en/stable/manual.html):
  polígonos planares, fronteiras, predicados, distâncias e validade.
- [NumPy array](https://numpy.org/doc/stable/reference/generated/numpy.array.html):
  representação dos dados numéricos.
- [Matplotlib Axes.set_aspect](https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.set_aspect.html):
  escala igual nos eixos das figuras.

As convenções de LE/TE e a política de folga acima são decisões do AeroOpt,
não definições aerodinâmicas atribuídas a essas bibliotecas.

## Geração NACA analítica e referências da v0.1.1

O gerador cria novas geometrias e nunca rediscretiza perfis importados.
Uma designação de quatro dígitos define `m=d1/100`, `p=d2/10` e
`t=(10*d3+d4)/100`. Para a família simétrica, exige-se prefixo `00`.
Na corda nominal unitária:

```text
h(x) = 5t [0.2969 sqrt(x) - 0.1260 x - 0.3516 x² + 0.2843 x³ - 0.1015 x⁴]
yc(x) = m/p² [2px-x²]                         se x <= p
yc(x) = m/(1-p)² [1-2p+2px-x²]                 se x > p
g(x) = dyc/dx; theta = atan(g)
P_upper = (x-h sin(theta), yc+h cos(theta))
P_lower = (x+h sin(theta), yc-h cos(theta))
```

Para `m=0`, `yc=g=0`. A espessura é aplicada normal à linha média, não
verticalmente. A formulação clássica tem espessura de fuga finita `2h(1)`;
para NACA 0012, `2h(1)=0.00252`. `closed_te=True` seleciona explicitamente o
coeficiente modificado `-0.1036`, obtido impondo `h(1)=0`. Essa modificação
é identificada nos resultados e não apresentada como a mesma geometria
clássica. O zero matemático da fuga modificada é representado exatamente.

Para N total ímpar, usam-se `(N+1)/2` estações de linha média de 0 a 1,
uniformes ou `x=(1-cos(beta))/2`, com beta uniforme em `[0,pi]`. A superfície
superior é percorrida ao contrário, seguida da inferior sem repetir a origem.
A referência nominal NACA `(0,0)` fica no índice `(N-1)/2`. As duas pontas de
fuga são mantidas; na fuga fechada são iguais. Não há transformação de
normalização oculta no gerador. O usuário escolhe explicitamente o índice
nominal ao construir `Airfoil`.

A origem nominal, o ponto de mínimo x da curva espessada e o ponto mais distante
de T são referências distintas em perfis cambados. No NACA 4412 com N=161
cossenoidal, o modo automático escolhe índice 79, e o nominal é 80. A diferença
é `0.00349061497 c`; a corda de referência automática difere em
`0.000300138689 c`. O mínimo x contínuo desse perfil está aproximadamente em
`(-0.000300247899, 0.003037524456)`. Nenhum deles substitui silenciosamente o
contrato nominal. Para S1223 publicado, usa-se explicitamente a amostra de
índice 45 em `(0.00005,0.00178)`: não é conhecida a LE exata da curva contínua.

As fontes das equações são as páginas PDAS de
[espessura](https://www.pdas.com/naca456thick4.html) e
[linha média](https://www.pdas.com/naca456mean2.html). A implementação é própria;
não se copiou um programa ou uma coleção de coordenadas externa para o pacote.
As verificações independentes incluem um cálculo Decimal, ordenadas publicadas,
integral exata de área simétrica e a referência escalar contínua descrita no
relatório `VALIDATION_V011.md`.
