# Método de painéis — AeroOpt v0.3.1

## Modelo e contrato

`PanelSolver` resolve escoamento potencial 2D, estacionário, incompressível e
invíscido ao redor de dois contornos fechados independentes. Há N fontes
constantes, uma por painel, e duas intensidades constantes de folha de vórtices,
uma por elemento. Todas as influências entre elementos são incluídas. Não há
camada limite, separação, transição, turbulência, estol, compressibilidade, asa
finita ou otimização. Reynolds permanece explícito no snapshot/hash, mas não
entra nas equações. Exige-se `mach=0` para selecionar inequivocamente o modelo;
Mach diferente de zero e extensões de condições são rejeitados, sem conversão.

O adaptador preserva `AerodynamicSolver.validate_input/prepare/execute/parse`,
`SolverRunner`, os esquemas v1 e a distinção SHA-256 científico/UUID operacional.
MSES e SU2 futuros não precisam adotar uma discretização ou campos de resultado
específicos deste método. Não há nova dependência de runtime.

```python
from pathlib import Path
from aeroopt.solvers import MomentReference, OperatingConditions, PanelSolver, SolverRunner

solver = PanelSolver(MomentReference(.25, 0, "clockwise", "q_inf*c_main^2"))
runner = SolverRunner(Path.cwd(), "outputs/panel_attempts")
result = runner.run(solver, geometry, OperatingConditions(6, 200000, 0),
                    min_clearance=.002)
```

`geometry` é uma montagem da API existente. A referência de momento é obrigatória:
posição global em cordas principais, sentido `clockwise` ou `counterclockwise`,
normalização exatamente `q_inf*c_main^2`. Não há referência física implícita.

## Painéis, orientação e singularidades

Cada segmento consecutivo da entrada vira um painel, sem redistribuição,
interpolação, suavização, remoção ou inversão de pontos. Uma fuga fechada com
extremidades exatamente iguais já contém todos os segmentos. Para fuga finita,
acrescenta-se somente o segmento de fechamento de `p_last` a `p_first`, a tampa
reta prevista em `SolverCase`. O snapshot/DAT original não muda. N inclui essas
tampas quando presentes. A colocação é no ponto médio de cada segmento.

Para tangente unitária `t=(tx,ty)`, sinal de área `s=+1` no percurso anti-horário
(CCW) e `s=-1` no horário (CW), a normal exterior é `n=s*(ty,-tx)`. A ordem de
entrada permanece preservada e sua orientação é registrada. O eixo y local do
painel aponta à esquerda de t: `l=(-ty,tx)`, portanto o exterior tem sinal
local `y=-s*|y|`. O comprimento L e o eixo x local seguem o segmento fornecido.

No ponto local `(x,y)`, defina `r1=hypot(x,y)`, `r2=hypot(x-L,y)`. A velocidade
de uma folha de fontes de intensidade unitária é

```text
u_s = [ln(r1) - ln(r2)] / (2*pi)
v_s = atan2(y*L, x*(x-L)+y*y) / (2*pi)
V_s = u_s*t + v_s*l
V_v = rotate_CCW_90(V_s)
```

O cálculo angular usa coordenadas escaladas antes dos produtos. A rotação de
90° define intensidade de vórtice **positiva CCW**, independentemente do percurso
do contorno. Os limites exteriores no ponto médio do próprio painel são
`V_s,self=n/2`, `V_v,self=s*t/2`; não se substitui y por um epsilon arbitrário.
As demais influências são integrais analíticas. Avaliar exatamente sobre outro
painel ou em sua extremidade é indefinido e gera rejeição. Fora do segmento,
um ponto colinear regular tem influência tangencial finita.

Desde a revisão v0.3.1, quando `|eta|<0.5`, o termo radial é calculado por
`atanh(eta)/(2*pi)`, com `eta=(r1²-r2²)/(r1²+r2²)`. O numerador é formado como
`L*(2*x-L)` nas coordenadas escaladas, sem subtrair raios quase iguais.
A identidade é exatamente a mesma integral; evita cancelamento de logaritmos
no campo distante. Para `|eta|>=0.5`, preservam-se os logaritmos para evitar
que eta arredondado a ±1 produza divergência falsa perto das extremidades.
Não se alteram pontos, sinais ou termos próprios.

Essas fórmulas podem ser obtidas integrando o campo de uma fonte pontual e
rotacionando-o para o vórtice. A formulação clássica está descrita em
[Mason, capítulo 4, seção 4.4](https://archive.aoe.vt.edu/mason/Mason_f/CAtxtChap4.pdf).
As convenções de orientação e sinais acima são as escolhas explícitas desta
implementação; não se copiou código de solver externo.

## Sistema acoplado e Kutta

As incógnitas são `z=(sigma_1,...,sigma_N,gamma_main,gamma_flap)`.
`V_inf=(cos(alpha),sin(alpha))`, com módulo unitário, conforme o contrato global
existente. Para cada colocação i:

```text
V_i = V_inf + sum_j sigma_j * V_s,ij
              + sum_e gamma_e * sum_(j in e) V_v,ij
V_i . n_i = 0
```

Cada linha de Kutta usa a primeira e a última colocação de **superfície** do
respectivo elemento, excluindo a tampa: `V_t,first + V_t,last = 0`, pois as
tangentes se opõem na fuga. A linha inclui fontes e vórtices dos dois elementos,
inclusive influências cruzadas e auto-influências. Existem N linhas de
impermeabilidade e duas de Kutta: sistema quadrado N+2, resolvido diretamente
por `numpy.linalg.solve`, sem mínimos quadrados, regularização ou pseudoinversa.

Para uma fuga aguda, a condição aproxima velocidades de saída iguais quando
os dois pontos médios se aproximam da fuga. Para uma fuga com espessura finita,
os dois pontos nunca coincidem: a mesma equação é uma **escolha de fechamento
aproximada**. A tampa recebe fontes e vórtices, impermeabilidade e integração
de pressão como qualquer outro painel. Não se impõe uma segunda condição nos
cantos da tampa, não se modela pressão de base nem se lança esteira explícita.

Alternativas seriam rejeitar toda fuga finita, alterar sua geometria para
fechá-la ou introduzir esteira/condição de base. Escolheu-se manter a tampa e
declarar a aproximação, preservando o contorno e o sistema N+2 solicitado.
Esse tratamento não foi validado experimentalmente, e seu refinamento pode
ser não monotônico. Uma cúspide de ângulo nulo também apresenta convergência
lenta, comprovada pelo caso exato de Joukowski. Kutta satisfeita em dois pontos
não certifica precisão da distribuição inteira nem do momento.

Na v0.3.1, qualquer elemento com tampa produz a flag legível por máquina
`limitations=["finite_te_closure_unverified"]` em `panel_output.json` e o
diagnóstico `FINITE_TE_UNVERIFIED` no resultado neutro. O parser exige a flag
correta para a geometria. A execução pode ter sucesso discreto, mas esses
coeficientes não devem ser usados para classificar projetos com fuga finita.
A revisão comparou subdivisões explícitas da mesma tampa e Kutta extrapolada;
nenhuma alternativa foi adotada automaticamente. Veja `VALIDATION_V031.md`.

## Pressão, forças, momento e circulação

No referencial adimensional `c_main=V_inf=1`, `Cp_i=1-V_t,i²`. A componente
normal é verificada antes de aplicar Bernoulli. Integração de pressão constante
em cada segmento usa `dF_i/(q_inf*c_main)=-Cp_i*n_i*L_i`. Como a força é
constante nesse painel reto, o momento é integrado exatamente nessa aproximação
com o braço do ponto médio.

```text
F = sum_i dF_i
CL = F . (-sin(alpha), cos(alpha))
CM_CCW = sum_i [(x_i-x_ref)*dFy_i - (y_i-y_ref)*dFx_i]
CM_CW = -CM_CCW
Gamma_e/(V_inf*c_main) = gamma_e * sum_(j in e) L_j
CL_circulation = -2 * sum_e Gamma_e/(V_inf*c_main)
```

CL usa a corda principal, sem dividir pela soma das cordas ou pela corda do
flap. CM usa `q_inf*c_main²`, mesmo quando a referência se desloca. Com nariz
a montante (x menor), momento positivo horário é nariz para cima; o sentido
CCW é o oposto. Não se escolhe automaticamente um desses sentidos.

CL publicado vem da pressão. Sua diferença em relação a Kutta–Joukowski é um
diagnóstico de consistência/refinamento, não eliminada substituindo CL pela
circulação. A componente de força na direção da corrente é registrada como
`streamwise_pressure_force_diagnostic`; **não é CD físico**. Um erro de
discretização ou uma aproximação de fuga não fornece arrasto viscoso.
`AerodynamicCoefficients.cd` e `panel_output.json.cd` permanecem `None`/`null`.

## Aceitação numérica e rejeições

A validação geométrica comum e a folga explícita solicitada pelo usuário
precedem a preparação. A política adicional do painel é fixa e entra em
`SolverIdentity.settings`, junto da versão do método, NumPy e referência CM:

| Verificação | Política da v0.3 |
|---|---|
| Comprimento de cada painel / corda própria | Estritamente maior que `1e-10` |
| Folga global | Pelo menos `1e-6 c_main`, além da restrição explícita do runner |
| Espessura de fuga / corda própria | No máximo 0.05 |
| Ordenação | LE interna, pelo menos dois segmentos por superfície; tangentes terminais voltadas para montante/jusante (`t_first.axis < -0.25`, `t_last.axis > 0.25`) e opostas (`dot < -0.5`) |
| Quantidade | No máximo 2000 painéis no total |
| Condicionamento | `cond_2(A) <= 1e10`, estimado por SVD; finito |
| Impermeabilidade | `max(abs(V_n))/V_inf <= 1e-9` |
| Kutta | `max_e abs(V_t,first+V_t,last)/V_inf <= 1e-9` |
| Erro retroativo | `||Az-b||inf/(||A||inf*||z||inf+||b||inf) <= 1e-9` |

Os limites de espessura, direção terminal, tamanho e folga são restrições
conservadoras de domínio/recursos, não fronteiras físicas validadas. A folga
não é confundida com gap aerodinâmico; registra-se também sua razão pelo maior
painel. Uma folga acima do limite pode continuar sub-resolvida e exige estudo
local de discretização. As tolerâncias geométricas preexistentes permanecem
intactas, com testes próprios. Painéis intrinsecamente pequenos são avaliados
na corda própria, sem confundir tamanho do flap com contato global.

Falha de condicionamento, singularidade, solução não finita ou resíduos
excessivos produz `accepted=false`, `not_converged`, sem Cp/CL/CM publicados.
Exceções inesperadas são `execution_failed`. Entradas rejeitadas são
`invalid_input`; saída inconsistente/truncada é `parsing_failed`. O parser
confere o caso, esquema, referência, dimensões, pressão versus velocidade,
circulação versus intensidade/comprimento e força/momento versus pressão.

`ConvergenceStatus.CONVERGED` significa aceitação do **problema discreto**.
Não certifica erro de discretização, convergência viscosa, validade experimental
ou precisão em dígitos. O condicionamento pode amplificar perturbações mesmo
com resíduos pequenos. O relatório de Joukowski mostra por que não se deve
usar esse estado isoladamente para avaliar CM próximo de zero.

## Execução e artefatos

A execução é uma biblioteca NumPy em processo, sem shell/processos externos.
Deadline é conferido antes/depois da montagem e após a solução/persistência,
mas não interrompe uma chamada BLAS bloqueante. Quem precisar de preempção
deve usar futuramente um worker por `AttemptContext.run_process`; o runner e
seu contrato de processos não foram alterados. Custo de armazenamento é O(N²),
de solução/SVD O(N³). Não há cache de resultado nem estado mutável no adaptador.

A montagem v0.3.1 projeta blocos de 64 painéis diretamente em coeficientes
normais/tangenciais. Não materializa os três tensores vetoriais N×N×2 da
estratégia anterior. As linhas normais são uma vista das primeiras linhas de A;
o custo dominante de armazenamento das duas matrizes é aproximadamente
`16*N²` bytes, além de blocos/temporários e workspace da álgebra linear.
A SVD completa de `cond_2` continua obrigatória em toda solução; limiares e
tratamento de falhas permanecem iguais. Versão do adaptador, estratégia de
montagem e tamanho do bloco fazem parte da identidade científica.

`panel_request.json` registra os segmentos exatos e sua relação com o snapshot;
`panel.log` registra o modelo. `panel_output.json` (esquema
`aeroopt.panel-output.v1`) inclui `case_id`, referência CM, aceitação/motivo,
diagnósticos e, apenas quando aceito, CL/CM/CD e os dois elementos em ordem
**main, flap**. Cada elemento contém centros, normais, comprimentos, Cp,
velocidades tangenciais/normais, fontes, intensidade de vórtice, circulação
CCW e fluxo líquido de fontes (diagnóstico de conservação). A última entrada
é a tampa, se `panel_count > surface_panel_count`. Os resíduos de Kutta também
seguem a ordem main/flap.

O resultado neutro guarda CL/CM/referência/CD ausente e o inventário SHA-256
desses artefatos. Leitura das distribuições usa `read_json` no arquivo da
tentativa. Artefatos de falhas e resultados parciais são preservados. Identidade
não contém tempo, caminho, PID ou UUID; a referência de momento e todos os
limites acima mudam o hash quando mudam a análise.

## Reprodução e limites científicos

```powershell
.\.venv\Scripts\python.exe examples/demo_panel_solver.py --output outputs/panel_demo_reproducao
.\.venv\Scripts\python.exe examples/verify_panel_solver.py --output outputs/panel_verification_reproducao
.\.venv\Scripts\python.exe examples/review_panel_solver.py --output outputs/panel_review_reproducao
.\.venv\Scripts\python.exe -m pytest -q
.\.venv\Scripts\python.exe -m pytest -q --reference-data=outputs/reference_data_v011
```

Escolha diretórios novos/vazios dentro do workspace. A demonstração utiliza os
DAT sintéticos existentes de 9/7 pontos, mantendo-os exatamente. O estudo cria
novas amostragens analíticas explícitas; não refina dados importados. Cilindro
e Joukowski exercitam o kernel de um corpo; o adaptador público permanece
restrito à montagem principal/flap. A referência Joukowski independente usa o
potencial exato do círculo e a transformação `w=z+1/z`, sem as influências ou
a solução linear de produção. Sua derivação e valores estão em
[VALIDATION_V030.md](VALIDATION_V030.md).

Não há comparação experimental viscosa nem certificação de desempenho. As
evidências apoiam consistência matemática em casos delimitados; não qualificam
geometrias arbitrárias, fugas grossas, folgas próximas do piso, alpha elevado
ou resultados para otimização. CL invíscido elevado com flap defletido não é
CL máximo nem prova de escoamento real aderido.

A [v0.3.2](VALIDATION_V032.md) preserva esta formulação e a identidade 0.3.1
do adaptador. Uma referência por imagens/conformação de dois corpos e estudos
HS/LV acrescentam evidência de refinamento; a reconstrução do aerofólio
publicado de Williams permanece não certificada. O módulo optativo
[`aeroopt.verification`](NUMERICAL_RELIABILITY.md) registra histórias e
classifica diferenças numéricas, sem alterar coeficientes, impor refinamento
automático ou fornecer limites rigorosos de erro. Fuga finita e momento
pequeno em cúspides continuam sujeitos às restrições da revisão.
