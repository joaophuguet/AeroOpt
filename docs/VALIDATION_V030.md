# Verificação do solver de painéis v0.3.0 — 2026-11-07

## Escopo, ambiente e reprodução

Implementação na branch `feature/panel-solver-v0.3`, a partir de
`ca3ad391375aa5b9d5969f7fc7a8e9171ad4ad0f`. O estado inicial de `main` estava
limpo. Não houve mudança de remoto, merge, push, rebase ou alteração dos dados
UIUC. O commit da entrega é informado na resposta final e no histórico Git.

Windows 11, CPython 3.13.15, NumPy 2.5.4. Código e metadados instalados do
AeroOpt estão em **0.3.0**. Todas as operações Python usaram exclusivamente
`.venv/Scripts/python.exe`; caches, artefatos e temporários ficaram no workspace.
Não foram adicionadas dependências de runtime nem instalações globais.

```powershell
.\.venv\Scripts\python.exe -m pytest -q
.\.venv\Scripts\python.exe -m pytest -q --reference-data=outputs/reference_data_v011
.\.venv\Scripts\python.exe examples/demo_panel_solver.py --output outputs/panel_demo_reproducao
.\.venv\Scripts\python.exe examples/verify_panel_solver.py --output outputs/panel_verification_reproducao
```

Os dois exemplos exigem saída nova/vazia no workspace. A versão final da
demonstração está em `outputs/panel_demo_v030_final/`; o estudo e a figura
revisada estão em `outputs/panel_verification_v030_reviewed/`. São artefatos
locais ignorados pelo Git, reprodutíveis pelos scripts versionados. As execuções
intermediárias foram preservadas. Nenhum dado aerodinâmico externo foi inventado
ou baixado para este estudo.

## Testes e auditoria

| Execução | Aprovados | Falhas | Ignorados | Warnings pytest | Tempo |
|---|---:|---:|---:|---:|---:|
| Base v0.2, antes da implementação | 278 | 0 | 21 | 0 | 30.71 s |
| v0.3, comando padrão completo | 337 | 0 | 21 | 0 | 66.15 s |
| v0.3, suíte completa com UIUC locais | 358 | 0 | 0 | 0 | 61.56 s |

São **59 testes novos** em `test_panel_method.py` e `test_panel_solver.py`.
Os 21 ignorados do modo padrão são exclusivamente os testes das coordenadas
UIUC fixadas por hash, que exigem `--reference-data`. Na execução com dados
locais, todos foram incluídos. Testes geométricos, contratos, falhas por etapa,
serialização, caminhos Windows, timeouts e encerramento de processos da v0.2
continuam passando. Nenhuma tolerância geométrica anterior foi alterada.

Cobertura nova: influências versus quadratura independente, limites singulares,
auto-influências nos dois percursos, sinal de circulação, cilindro, simetria,
inclinação de sustentação, Joukowski, Kutta independente, interação/distanciamento,
sensibilidade de posição/deflexão, refinamento, tampa finita, mudança de
referência/sentido de momento, identidade/repetições, preservação de entradas,
coeficientes ausentes em falhas, corrupção/truncamento de saída, deadline
cooperativo e CLIs reprodutíveis com proteção contra sobrescrita.

Foram reimportados os 13 resultados bem-sucedidos da demonstração e do estudo
(1 + 12 tentativas), conferindo os hashes SHA-256 e tamanhos de 130 artefatos.
As figuras de geometria/Cp e refinamento foram inspecionadas visualmente.
`pip check` informou ausência de dependências quebradas.

## Influências analíticas e cilindro

A integral de fonte em segmento é comparada com quadratura Gauss–Legendre de
128 pontos em três orientações. O campo deve permanecer igual ao inverter o
segmento. Testam-se os limites exteriores/interiores `(u,v)` em `(0.4,±1e-12)`
para segmento unitário: `u=ln(0.4/0.6)/(2*pi)`, `v→±0.5`. Os pontos exatamente
sobre a folha/extremidades são explicitamente rejeitados. Os termos próprios
exteriores são inseridos analiticamente e verificados para CW e CCW.

O cilindro não possui Kutta: o teste resolve apenas o bloco de fontes do mesmo
kernel. Para um polígono regular, compara-se Cp no ângulo de cada ponto médio
com `Cp(theta)=1-4*sin(theta)^2`. A simetria especial dessa discretização produz
acordo quase exato nos ângulos de colocação; não significa precisão exata para
um contorno arbitrário ou ao longo de todos os segmentos.

| Painéis do cilindro | Erro máximo de Cp | Resíduo normal máximo |
|---:|---:|---:|
| 24 | 1.78e-15 | 7.77e-16 |
| 48 | 4.00e-15 | 8.88e-16 |
| 96 | 5.77e-15 | 1.78e-15 |

## Simetria e limite de aerofólio fino

Dois NACA 0012 analíticos, fuga fechada, 161 pontos por elemento, razão 0.3,
LE do flap `(1.05,0)`, deflexão e alpha zero: `CL=3.11e-17`,
`CM_CW=6.51e-18`, circulações CCW `(-3.82e-17,4.70e-18)`. As duas condições
de Kutta têm resíduo zero nesta execução; o máximo normal é `8.88e-16`.
Casos de um elemento com fuga aguda e finita também verificam Cp espelhado.

O limite fino usa NACA 0002, fuga fechada, incidências ±2°, inclinação secante
por radiano. `2*pi=6.283185307...` é o limite de teoria fina, não uma solução
exata de um perfil com espessura de 2%.

| Painéis | dCL/dalpha, rad⁻¹ | Diferença relativa de 2*pi |
|---:|---:|---:|
| 80 | 6.48602787 | +3.228% |
| 160 | 6.41098771 | +2.034% |
| 320 | 6.39369449 | +1.759% |

O sinal da sustentação e sua tendência são consistentes. A diferença remanescente
inclui espessura finita e discretização; não é eliminada ajustando tolerâncias.

## Referência exata independente de Joukowski

`examples/_panel_reference.py` implementa o potencial analítico de um círculo,
sem chamar a integração de influências ou o sistema linear de produção.
O círculo tem centro `z0=-0.08`, raio `R=1.08`, passando pelo ponto crítico
`z=1`. A transformação é `w=z+1/z`; a corda mapeada é
`C=2+(1+2*0.08)+1/(1+2*0.08)`. Normalizam-se as coordenadas por C e desloca-se
a LE para zero. A malha é gerada explicitamente em theta uniforme no círculo,
mantendo a cúspide repetida no fim. Não há arquivo importado rediscretizado.

Para `alpha=4°`, o potencial complexo e sua derivada dão

```text
Gamma = -4*pi*R*sin(alpha)                 # CCW positivo
dW/dz = exp(-i*alpha) - R²*exp(i*alpha)/(z-z0)²
         + Gamma/[2*pi*i*(z-z0)]
dW/dw = (dW/dz)/(1-1/z²)
Cp = 1-|dW/dw|²
CL = -2*Gamma/C
```

O momento de referência é a integral da pressão exata, positivo CCW em
`(0.25,0)`. Quadraturas periódicas de ponto médio com 16384 e 32768 amostras
concordam em CM dentro de `1e-12`; a integração de pressão reproduz CL exato
também dentro de `1e-12`. Isso verifica a referência numérica independente,
sem estabelecer um limite rigoroso de erro.

Referências: `CL=0.47075898801417243`,
`Gamma/(V_inf*c_main)=-0.23537949400708622`,
`CM_CCW=+0.001240656232512379`.

| Painéis | CL por pressão | Erro absoluto CL | CM CCW | Erro absoluto CM |
|---:|---:|---:|---:|---:|
| 80 | 0.41275967 | 0.05799932 | -0.02295916 | 0.02419982 |
| 160 | 0.43330686 | 0.03745212 | -0.01527661 | 0.01651726 |
| 320 | 0.44795665 | 0.02280234 | -0.00921545 | 0.01045611 |
| 640 | 0.45740361 | 0.01335538 | -0.00505882 | 0.00629948 |
| 1920 | 0.46531652 | 0.00544247 | -0.00140820 | 0.00264886 |

O erro de CL na última malha é **1.156%**; o de circulação é **0.424%**.
Os erros de CL, circulação e CM diminuem em todas as malhas. Entretanto,
**o sinal do pequeno momento ainda é incorreto em 1920 painéis**. Não se
declara CM preciso nessa cúspide. `cond_2` cresce de aproximadamente 854 para
95792; os resíduos normais continuam <=`2.33e-13` e Kutta <=`3.42e-14`.
Resíduos baixos não detectam esse erro de discretização.

O teste de momento inicialmente não alcançou o alvo absoluto `0.003` com 640
painéis (`0.00629948` observado). O caso foi reproduzido e refinado a 1920,
mantendo o alvo; não se enfraqueceu esse critério para aprovar a malha original.
A evidência sustenta convergência lenta de painéis constantes nessa cúspide,
e não adequação do método para prever momentos pequenos sem estudo próprio.

## Principal + flap: interação e refinamento

NACA 0012 principal e NACA 4412 flap, ambos analíticos de fuga fechada;
referências nominais LE explícitas, razão 0.3, LE do flap `(1.02,-0.04)`,
deflexão horária 25°, alpha 6°, Mach 0 e Reynolds registrado 200000.
CM é positivo horário em `(0.25,0)` e normalizado por `q_inf*c_main²`.
As amostragens cossenoidais são geradas explicitamente para cada malha.

| Painéis por elemento | Total N | CL pressão | CM horário | \|CL pressão − CL circulação\| |
|---:|---:|---:|---:|---:|
| 40 | 80 | 3.12020216 | -0.77907231 | 0.05446271 |
| 80 | 160 | 3.14061203 | -0.79406030 | 0.02949275 |
| 160 | 320 | 3.14942923 | -0.80195983 | 0.01525084 |
| 320 | 640 | 3.15316882 | -0.80591799 | 0.00774253 |

Na última malha, `Gamma_main=-1.23451219`, `Gamma_flap=-0.34594348`,
`CL_circulation=3.16091135`, `cond_2=231.05`; os dois resíduos de Kutta são
`(-8.88e-16,3.997e-15)` e a impermeabilidade máxima é `1.92e-14`.
As diferenças entre malhas consecutivas diminuem para CL e CM; a última
diferença é 0.00373959 em CL e 0.00395815 em CM. O desacordo pressão/circulação
também diminui e ainda é aproximadamente 0.246% de CL. Não há referência exata
multielemento neste estudo: consistência/refinamento não certificam erro zero.

Os testes comparam cada circulação acoplada à solução de seu elemento isolado:
ambas mudam mais de 0.05 na configuração de teste. Ao mover o flap a `(100,0)`,
as circulações retornam às isoladas dentro de 1%. A independência das duas
intensidades é verificada com dois valores distintos e ambas as equações de
Kutta satisfeitas, sem impor uma diferença física arbitrária entre intensidades.

Sensibilidade com 160 painéis por elemento e os demais dados acima:

| Configuração | CL | CM horário |
|---|---:|---:|
| Flap a 0°, LE `(1.02,-0.04)` | 1.28013531 | -0.19794033 |
| Flap a 20°, mesma LE | 2.80285220 | -0.69012972 |
| Flap a 20°, LE `(1.20,-0.10)` | 2.46389216 | -0.88566282 |

Esses valores invíscidos não preveem adesão real ou estol. Não se interpreta
CL elevado como CL máximo, nem se recomenda uma configuração de projeto.

## Fuga finita: fechamento aproximado

Repetiu-se a montagem de 25° com NACA clássicos de fuga finita, acrescentando
uma tampa por elemento. Kutta é imposta nas duas superfícies adjacentes,
excluindo a tampa; impermeabilidade e pressão incluem a tampa. Esse modelo
não fornece pressão de base ou separação nos cantos.

| Total N, incluindo 2 tampas | CL pressão | CM horário | \|CL pressão − CL circulação\| |
|---:|---:|---:|---:|
| 82 | 3.01715919 | -0.76006841 | 0.03922681 |
| 162 | 3.14220344 | -0.80305378 | 0.01312442 |
| 322 | 3.14274269 | -0.80470144 | 0.00608997 |
| 642 | 3.13160224 | -0.80225805 | 0.00266456 |

O desacordo pressão/circulação diminui, mas CL e CM **não convergem
monotonicamente** nessa sequência. A variação de CL de 322 para 642 painéis é
0.01114045, maior que a variação anterior. Não se usa a proximidade dos dois
valores intermediários como certificado de convergência. A componente de força
na corrente na última malha é `-0.00220262`, somente diagnóstico; CD segue null.
O tratamento de fuga finita precisa de referência independente e estudo
adicional antes de classificação de candidatos.

## Exemplo com a geometria de demonstração existente

Os DAT de 9 e 7 pontos foram importados sem alterar amostras. Razão 0.3,
LE do flap `(1.10,-0.05)`, deflexão 25°, alpha 6°, Re registrado 200000,
Mach 0, folga requerida 0.002. Com as tampas, são 16 painéis e 18 incógnitas.

Resultado de execução: `CL=2.1977484228`, `CM_CW=-0.7256996941`, `CD=null`.
Esses polígonos muito grosseiros demonstram a integração com o runner e a
publicação de Cp; não são um estudo convergido nem perfis com desempenho
experimental conhecido. O gráfico exibe marcadores desconectados nos pontos
médios, incluindo a tampa, sem sugerir uma interpolação da entrada.

## Falhas exercitadas e limites restantes

Geometria inválida, ordenação incompatível, painel nulo/sub-resolvido,
fuga fora do domínio, Mach não suportado e folga abaixo de `1e-6 c_main`
são rejeitados sem modificar pontos. Sistemas singulares, `cond_2>1e10`,
erro de álgebra linear, solução não finita e resíduos excessivos não publicam
coeficientes. Corrupção de caso, Cp, malha, circulação, CM ou tentativa de
publicar CD produz falha de parsing. Falhas numéricas previsíveis preservam
diagnósticos e `not_converged`; exceções preservam traceback e evidências.

A folga numérica mínima é uma restrição conservadora, sem validação de precisão
perto desse piso. O limite de painéis é 2000 no total, com memória O(N²),
SVD/solução O(N³). Não há remalhamento automático nem estudo universal de
qualidade. Timeout da biblioteca é cooperativo; BLAS não pode ser preemptido.

A primeira atualização de metadados por pip falhou porque setuptools não estava
instalado na `.venv`, e a tentativa com build isolado não tinha rede disponível.
Resolveu-se offline copiando para `.tmp/build_wheels_v030/` o wheel setuptools
84.2.1 já existente no cache local e usando `--no-index --find-links` no build
isolado. Não se instalou dependência global nem se recriou a `.venv`.
Uma primeira execução do relatório também revelou um escalar NumPy incompatível
com o JSON estrito; a inclinação passou a ser convertida explicitamente para
float Python, coberta agora pelo teste completo da CLI.

O conjunto fornece evidência de consistência matemática e comportamento de
falhas do solver potencial discreto. Não constitui validação experimental,
garantia de precisão para todas as geometrias ou prontidão para otimização
aerodinâmica. Momento próximo de zero, cúspides e fugas finitas permanecem
limitações quantitativamente demonstradas, além da ausência deliberada de
física viscosa, arrasto e estol.
