# Revisão científica e numérica v0.3.1 — 2026-11-07

## Parecer técnico

**O solver ainda não está qualificado para otimização aerodinâmica automática,
mesmo preliminar.** É utilizável em estudos exploratórios invíscidos controlados,
com referência e refinamento próprios. A revisão confirma a formulação e as
convenções nos casos exercitados, corrige um defeito de precisão de influência
e reduz o custo de montagem/memória. Não resolve a precisão de CM pequeno na
cúspide de Joukowski nem estabelece um fechamento confiável de fuga finita.
Não há validação experimental, previsão de arrasto viscoso ou estol.

| Questão | Constatação | Decisão |
|---|---|---|
| Sinais, normais e momento | Derivação independente, quadratura de influências e Blasius concordam | Preservar convenções |
| Campo distante | Subtração de logaritmos perde precisão relativa | Corrigir com identidade equivalente; regressão Decimal |
| CM de Joukowski | Erro dominante na velocidade discreta; não é erro de sinal ou referência | Preservar método e declarar limitação de precisão |
| Fuga finita | Sensível à tampa e à posição das colocações, com condicionamento baixo | Manter aproximação com flag obrigatória de não verificada |
| Interação multielemento | Diferenças entre Hess–Smith e vorticidade linear diminuem com refinamento | Verificação numérica adicional, sem alegar referência exata |
| Custo | Menos tensores; SVD domina em malhas grandes | Montagem em blocos; manter SVD e limites de aceitação |

## Branch, ambiente e reprodução

Base: `0b8118f17b8c63bea588556677491b19222f712a`, árvore inicialmente limpa,
branch `feature/panel-solver-v0.3`. Não houve merge em `main`, push ou mudança
de remoto. O commit final da revisão é informado na resposta e no histórico.

Windows 11, Python 3.13.15, NumPy 2.5.4. Código e metadados do AeroOpt: 0.3.1.
Usou-se exclusivamente `.venv/Scripts/python.exe`, sem recriar o ambiente ou
adicionar dependências de runtime. A atualização editável usou `--no-deps
--no-index --find-links .tmp/build_wheels_v030`, com o wheel de build já em
cache. `pip check` não encontrou requisitos quebrados.

```powershell
.\.venv\Scripts\python.exe -m pytest -q
.\.venv\Scripts\python.exe -m pytest -q --reference-data=outputs/reference_data_v011
.\.venv\Scripts\python.exe examples/review_panel_solver.py --output outputs/panel_review_reproducao
```

Escolha uma pasta nova/vazia no workspace. O estudo completo final fica em
`outputs/panel_review_v031_release/`; execuções exploratórias anteriores foram
preservadas. `--quick` serve apenas ao teste da CLI, sem substituir este estudo.
`review.json` guarda métricas, ambiente, configurações e hashes do código;
`cases/` guarda coordenadas, alpha, referência de momento, proveniência e saídas
completas. `sources/` preserva os seis arquivos de implementação do estudo;
`artifacts.json` inventaria SHA-256/tamanho dos arquivos, excluindo a si mesmo.
`input_sha256` dos casos de pesquisa identifica a entrada canônica, não substitui
o `SolverCase.case_id` do runner. Tempos não entram nesses hashes.

São geometrias analíticas geradas explicitamente. Nenhum DAT importado ou
original UIUC foi reordenado, interpolado ou alterado. Não foram incorporadas
coordenadas externas ao Git. Saídas e cópias de código geradas ficam em outputs.

## Auditoria independente da formulação

Para uma fonte pontual unitária, `dV=(r/|r|²)ds/(2*pi)`. Num painel de
comprimento L com x ao longo da tangente e y à esquerda, integrar de 0 a L dá

```text
u_s = integral (x-s)/[(x-s)^2+y^2] ds/(2*pi)
    = ln(r1/r2)/(2*pi)
v_s = integral y/[(x-s)^2+y^2] ds/(2*pi)
    = atan2(y*L, x*(x-L)+y*y)/(2*pi)
V_v = (-V_s,y, V_s,x)                  # vórtice CCW positivo
```

Para percurso de sinal de área `s=±1`, `n=s*(ty,-tx)` é exterior. O limite
normal de fonte exterior é +1/2, logo `V_s,self=n/2` e
`V_v,self=s*t/2`. Inverter o percurso inverte t e s, preservando o campo físico
e Gamma CCW. A revisão verifica os dois percursos e ambos os tipos de fuga.
Quadraturas independentes integram os campos pontuais, não os coeficientes do
kernel. Não foi encontrado erro de sinal nessas relações.

O sistema contém N intensidades de fonte e uma intensidade de vórtice por
elemento. Cada coluna de vórtice soma somente os painéis do elemento associado,
mas afeta **todas** as colocações. Cada linha Kutta soma as velocidades
tangenciais das duas colocações terminais de superfície do seu elemento;
inclui todos os campos cruzados. São N+2 equações para dois elementos, sem
forçar intensidades iguais. A tampa não é usada como colocação Kutta.

De `dF/(q*c_main)=-Cp*n*ds`, seguem `CL=F.(-sin(alpha),cos(alpha))` e
`CM_CCW=sum((x-x_ref)*dFy-(y-y_ref)*dFx)`. Para deslocamento da referência
`(dx,dy)`, `CM_new=CM_old-dx*Fy+dy*Fx`; CM horário é seu negativo. A integração
do momento no ponto médio é exata para **pressão constante no segmento reto**;
não é exata para o campo contínuo sobre o perfil curvo.

A folha de vórtices tem `Gamma_e=gamma_e*perimetro_e`; sua expansão distante
dá `CL_Gamma=-2*sum(Gamma_e)/(V_inf*c_main)`. O fluxo líquido das fontes é
diagnóstico: a impermeabilidade em colocações não garante impermeabilidade
contínua nem fluxo integrado exatamente nulo. Portanto CL por pressão e por
circulação podem diferir em malha finita. Ambos são preservados e comparados;
não se substitui um pelo outro para ocultar o desacordo. O arrasto físico
continua indisponível, `cd=None`/null.

## Correção confirmada: cancelamento no campo distante

A expressão anterior `log(r1)-log(r2)` subtrai números quase iguais.
Usa-se agora, em coordenadas escaladas,

```text
eta = L*(2*x-L)/(r1^2+r2^2)
log(r1/r2) = atanh(eta)
```

Para `|eta|<0.5`, a forma atanh evita cancelamento. Perto das extremidades,
mantêm-se os logaritmos para evitar eta arredondado a ±1. Isso não altera a
integral, os termos próprios ou a geometria. A referência usa Decimal com 60
dígitos para um painel unitário horizontal observado em `(R,0)`.

| R | Erro relativo anterior | Erro relativo após correção |
|---:|---:|---:|
| 1e6 | 1.84e-10 | 3.33e-16 |
| 1e9 | 8.22e-8 | <=2.22e-16 |
| 1e12 | 1.69e-3 | <=2.22e-16 |

A regressão com tolerância relativa `1e-13` falharia na implementação anterior.
Também se exercitam pontos regulares muito próximos das extremidades. Essa
correção melhora a robustez da influência; **não corrige o CM de Joukowski**:
nos casos usuais as alterações em CL/CM são da ordem do arredondamento.

## Joukowski: confirmação da referência e decomposição do erro

Mantém-se a geometria anterior: círculo `z0=-0.08`, raio `R=1.08`, transformação
`w=z+1/z`, alpha 4°, corda `C=2+1.16+1/1.16`, referência `(0.25,0)`, CM CCW.
O teste de pressão analítica foi confirmado por **Blasius em forma fechada**.
Com `U=exp(-i*alpha)` e circulação dimensional `G=-4*pi*R*sin(alpha)`,

```text
dW/dw = U + A/w + B/w^2 + ...
A = G/(2*pi*i)
B = U - R^2*conj(U) + i*G*epsilon/(2*pi)
M_origin/rho = pi*imag(A^2+2*U*B)
CM_ref = 2*[pi*imag(A^2+2*U*B) + G*(x_ref*cos(alpha)+y_ref*sin(alpha))]/C^2
```

Aqui x_ref e y_ref estão no plano w antes de normalizar a corda. Essa fórmula
não integra Cp nem chama um solver. Concorda com a quadratura independente
de pressão dentro de `2e-13`, inclusive para incidências negativas e zero:
`CL=0.4707589880141724`, `CM_CCW=+0.00124065623251235`.

| Painéis | CM Hess–Smith | CM com Cp analítico sobre a mesma poligonal | CM de vórtices lineares |
|---:|---:|---:|---:|
| 80 | -0.02295916 | +0.00131896 | +0.00225248 |
| 160 | -0.01527661 | +0.00126027 | +0.00171949 |
| 320 | -0.00921545 | +0.00124556 | +0.00147309 |
| 640 | -0.00505882 | +0.00124188 | +0.00135509 |
| 1920 | -0.00140820 | +0.00124079 | Rejeitado por condicionamento |

**Confirmado:** erro de sinal/referência não explica a discrepância; o mesmo
polígono, normais e integração recuperam o sinal usando pressão analítica.
Em 160 painéis, o erro geométrico/quadratura desse diagnóstico é `1.96e-5`,
contra erro de CM Hess–Smith `1.65e-2`. O erro dominante está no campo de
velocidades/pressões da solução discreta, não na multiplicação do braço.

**Ensaios diagnósticos, não correções adotadas:** em 320 painéis, impor a
circulação exata e resolver as equações normais dá CM aproximadamente
`-0.00633984`, ainda com sinal errado, violando a Kutta terminal original em
cerca de 0.261. Portanto o erro de circulação não explica sozinho o momento.
Extrapolar linearmente as velocidades dos dois painéis terminais de cada lado
até a fuga dá CM aproximadamente `-0.00942940`, ligeiramente pior que o atual.
Esses valores são marcados como diagnósticos e não aceitos pelo adaptador.

Também se gera explicitamente `theta=beta-0.5*sin(beta)`, concentrando mais
amostras na cúspide. Em 1920 painéis, CM melhora para cerca de `-0.00054018`,
mas continua com sinal errado; o erro absoluto ainda é `0.00178084`.
Não se redistribuiu nenhum arquivo de entrada nem se adotou esse viés em
produção. As tabelas completas incluem CL, Gamma, condicionamento e resíduos.

**Interpretação sustentada:** discretização constante, geometria quase
coincidente junto à cúspide e condição Kutta discreta participam do erro.
A troca da distribuição e da formulação muda fortemente a convergência.
Não foi estabelecida uma decomposição causal única entre essas três fontes.
Não se encontrou justificativa quantitativa para alterar a Kutta de produção.

**Limitação explícita:** no limite prático de 2000 painéis, não se consegue
qualificar o sinal de CM da ordem de `1e-3` nessa cúspide. O erro absoluto
medido da malha uniforme de 1920 painéis é `0.00264886` e o de CL é 1.156%.
Há convergência em direção à referência, mas isso não constitui precisão
suficiente de momento pequeno. Não se propõe `0.00265` como limite universal
de erro para outros perfis. Nenhuma tolerância anterior foi relaxada.

## Referência independente multielemento

A busca priorizou fontes primárias:

- [Williams, ARC R&M 3717, 1971/1973](https://reports.aerade.cranfield.ac.uk/bitstream/handle/1826.2/2993/arc-rm-3717.pdf?isAllowed=y&sequence=1):
  caso exato de dois aerofólios por imagens e duas transformações Kármán–Trefftz.
  O relatório foi localizado, com equações e tabelas impressas. Não se reconstruiu
  a dupla transformação/seleção de ramos nem se digitalizaram coordenadas para
  tratá-las como uma referência de alta resolução. Essa comparação permanece
  pendente; não se afirma que não exista benchmark publicado adequado.
- [Suddhoo e Hall, 1985](https://doi.org/10.1017/S0001924000096767): casos exatos
  de três/quatro elementos. O resumo foi acessível, sem conjunto numérico
  completo utilizável nesta execução.
- [Cox, 2011](https://digitalcommons.calpoly.edu/aerosp/28/) e a
  [descrição do método de vórtices lineares da Virginia Tech](https://www.engapplets.vt.edu/fluids/vpm/vpminfo.html)
  apoiam a escolha da alternativa de formulação. Não foram copiados programas,
  inferidos dados ausentes ou usados coeficientes publicados sem geometria.

Adotou-se a **segunda opção solicitada**: implementação própria e separada de
vorticidade linear contínua por segmentos, com N+1 valores nodais por contorno,
nenhuma fonte e nenhuma intensidade constante global por elemento. Cada fuga
tem `gamma_first+gamma_last=0`. Não importa `aeroopt.solvers`, não usa suas
influências ou montagem, e continua funcionando em teste quando essas funções
de produção são substituídas por exceções.

Sua integral complexa parte diretamente do vórtice pontual CCW. Com
`q=(z-z_start)/t`, `I0=log(q/(q-L))/t` e
`I1=[(q/L)*log(q/(q-L))-1]/t`, as duas funções de forma têm velocidades
`(I0-I1)/(2*pi*i)` e `I1/(2*pi*i)`. Os coeficientes são comparados com
quadratura de vórtices pontuais. Limites próprios CW/CCW são derivados
separadamente. Séries no campo distante evitam cancelamento em ambos os
momentos da integral. Condicionamento e resíduos também são verificados.

Essa referência tem limitações próprias: para Joukowski uniforme a condição
cresce para aproximadamente `6.8e9` em 640 painéis e `3.2e11` em 1920. O último
caso é rejeitado por exceder `1e10`, sem publicar coeficientes de referência.
Cp nas colocações mais próximas da cúspide não está qualificado, mesmo onde
os momentos integrados melhoram. Por isso não se usa esse método como verdade
exata nem como certificado geral de Cp. Os casos NACA de ângulo de fuga finito
e espessura de fuga **zero** têm condicionamento bem menor.

### Comparação controlada de dois elementos

Geometria analítica da v0.3: NACA 0012/4412, fuga fechada, razão 0.3, LE do
flap `(1.02,-0.04)`, deflexão 25°, alpha 6°. Referências LE nominais explícitas.
Ambos os métodos recebem exatamente os mesmos arrays globais e comparam Cp
nas mesmas colocações, sem interpolação. CM abaixo é CCW em `(0.25,0)`.

| Total de painéis | CL HS − CL linear | CM HS − CM linear | RMS Cp principal | RMS Cp flap |
|---:|---:|---:|---:|---:|
| 160 | +0.02066858 | -0.01255257 | 0.228478 | 0.132815 |
| 320 | +0.01087220 | -0.00614771 | 0.114610 | 0.061232 |
| 640 | +0.00568148 | -0.00300475 | 0.058252 | 0.028948 |
| 1280 | +0.00292681 | -0.00147815 | 0.030060 | 0.013999 |

Todas essas diferenças diminuem. Em 1280 painéis, CL HS é aproximadamente
3.15478454 e CL linear 3.15185772; CM HS 0.80787128 e CM linear 0.80934943.
O condicionamento da referência linear é aproximadamente `3.52e5`, abaixo
do limite, e suas circulações são aproximadamente `(-1.232300,-0.345764)`.
As circulações HS são `(-1.233499,-0.345843)`, com diferenças absolutas
`(0.00119916,0.00007885)`. O RMS Cp ponderado por comprimento na mesma malha
é `(0.015498,0.007885)` para principal/flap.
Os JSON preservam as circulações de cada método, Cp completo e também RMS
ponderado por comprimento. Acordo entre dois métodos não fornece erro absoluto
garantido: compartilham o mesmo modelo invíscido e a mesma aproximação poligonal.

Sensibilidade em 640 painéis totais; as posições são LE global do flap:

| LE do flap / deflexão | CL HS | CL linear | CM HS CCW | CM linear CCW |
|---|---:|---:|---:|---:|
| `(1.02,-0.04)` / 25° | 3.153169 | 3.147487 | 0.805918 | 0.808923 |
| `(0.98,-0.04)` / 25° | 3.127112 | 3.124246 | 0.745964 | 0.750507 |
| `(1.10,-0.04)` / 25° | 2.981090 | 2.974358 | 0.906920 | 0.908306 |
| `(1.02,-0.07)` / 25° | 3.021262 | 3.016122 | 0.801699 | 0.804383 |
| `(1.02,-0.04)` / 15° | 2.443951 | 2.440088 | 0.576081 | 0.578810 |

Isso exercita afastamento vertical, deslocamento longitudinal/sobreposição
nominal e deflexão. O deslocamento longitudinal é `x_LE_flap-1`; seu negativo
é a sobreposição nominal nesta definição. A distância mínima global é guardada
separadamente, sem chamá-la de gap aerodinâmico padronizado. Não há validação
de uma faixa universal de gap/overlap, de incidência alta ou de estol.

## Fuga finita: estudo e recomendação

Mantém-se a montagem anterior, com NACA clássicos de fuga finita. A tampa
reta representa a mesma região sólida quando subdividida em 1, 4 ou 8 painéis;
os segmentos de superfície e as colocações Kutta não mudam entre essas três
variantes na mesma linha. A variação abaixo é numérica, sem alterar o perfil.

| Painéis de superfície por elemento | CL, 1 painel de tampa | CL, 4 painéis | CL, 8 painéis |
|---:|---:|---:|---:|
| 40 | 3.017159 | 3.027461 | 3.028641 |
| 80 | 3.142203 | 3.164179 | 3.165972 |
| 160 | 3.142743 | 3.167403 | 3.170059 |
| 320 | 3.131602 | 3.161500 | 3.164805 |
| 640 | 3.118663 | 3.158319 | 3.162336 |

Na última linha, a troca de 1 para 8 painéis altera CL em **0.04367286** e
CM CCW de 0.79898282 para 0.81126205. As condições são aproximadamente 320 e
317, muito distantes de `1e10`, com resíduos pequenos. Não é evidência de
instabilidade da solução linear; é sensibilidade da discretização/fechamento.
A diferença entre 4 e 8 tampas ainda é 0.00401667 em CL. Não há referência
exata de pressão de base para decidir qual curva representa melhor o físico.

Também se geram novas amostras NACA com `x=1-(1-x_cosine)^p`. Com 320 painéis
de superfície por elemento e uma tampa, `p=0.75,1,1.25` produz CL de
`3.126995,3.131602,3.112193`. A pressão e a Kutta terminal são sensíveis à
resolução local perto dos cantos; refinar somente a superfície não refina a
representação de pressão sobre a tampa.

A alternativa de extrapolar a velocidade das duas colocações de cada lado
até os cantos também foi comparada. Por exemplo, em 1280 painéis de superfície
e uma tampa por corpo, CL passa de 3.118663 para aproximadamente 3.115814;
subdividir a tampa em oito painéis produz aproximadamente 3.162379 com essa
alternativa. Ela não elimina a dependência da tampa. Seus resíduos de equação
e sua diferença da Kutta original ficam registrados como diagnósticos; não é
uma solução aceita pelo adaptador atual.

**Recomendação:** não adotar ajustes para tornar as curvas monotônicas.
Preservar a aproximação existente com `finite_te_closure_unverified` obrigatório
no artefato e diagnóstico `FINITE_TE_UNVERIFIED` no `SolverResult`. Geometria
válida e solução discreta bem condicionada não qualificam esses coeficientes
para classificar configurações de fuga finita. O modelo não impõe condições
independentes nos dois cantos, não calcula pressão de base e não modela esteira
ou separação. A raiz quantitativa completa da convergência ainda está aberta.

## Custo e melhoria implementada

O benchmark compara a estratégia vetorial v0.3.0 com a nova projeção em blocos,
usando os mesmos coeficientes corrigidos, geometria e backend. A estratégia
anterior preservada no harness **não** é usada como referência independente.
São três repetições por caso, mediana, após aquecimento básico da álgebra
linear. Medições finais são sequenciais, sem pytest concorrente. Timings excluem
runner, validação geométrica, parsing e I/O; não são promessas de throughput.

| N total | Núcleo anterior, s | Núcleo em blocos, s | Pico rastreado anterior, MiB | Pico em blocos, MiB |
|---:|---:|---:|---:|---:|
| 80 | 0.0176 | 0.00280 | 0.453 | 0.683 |
| 160 | 0.0552 | 0.0283 | 1.783 | 1.961 |
| 320 | 0.1418 | 0.0926 | 7.080 | 4.930 |
| 640 | 0.3561 | 0.2701 | 28.220 | 12.977 |
| 1280 | 1.2797 | 1.1560 | 112.688 | 38.446 |
| 1920 | 3.8199 | 3.2845 | 253.405 | 85.667 |

Os números acima são do lote final `release/`, que preserva todas as amostras
brutas. Um lote isolado anterior (`outputs/panel_review_v031_final/`) mediu
aproximadamente 3.18 s em ambas as montagens com N=1920. Essa variação impede
prometer um ganho fixo de throughput. Mediu-se memória com tracemalloc em passagem **separada** das
medições de tempo, evitando favorecer blocos pelo overhead de tracing sobre o
loop Python antigo. Esse pico é de alocações rastreadas por chamada, não RSS
total: workspace nativo não rastreado, interpretador e caches BLAS ficam fora.
Em N pequeno os blocos têm overhead maior; em 1920 a redução é cerca de 66%.

A memória persistente dominante passa a duas matrizes escalares, aproximadamente
`16*N²` bytes, mais O(N*64) de blocos e workspace. A montagem ainda é O(N²);
solução direta e SVD continuam O(N³). Em 1920 a SVD levou cerca de 2.25 s e
continua dominando o custo; o ganho total variou de praticamente zero a cerca
de 14% entre os dois lotes. Em 640 houve redução de aproximadamente 24% no
núcleo do lote final. Mil geometrias com esses tempos custariam cerca de
4.5 min em N=640, 19.3 min em N=1280 e 55 min em
N=1920, **antes** dos custos excluídos. Isso não inclui estudos de refinamento
por geometria nem qualquer otimização implementada.

A maior diferença entre matrizes foi aproximadamente `7.11e-15`; diferenças
de CL/CM no benchmark foram <=`9e-16`. Regressões com CW/CCW e tampas comparam
as equações, soluções e coeficientes. A versão/estratégia/tamanho de bloco entram
na identidade; os inputs globais permanecem imutáveis. Os gates originais de
condicionamento, singularidade, finitude, resíduos e falhas foram mantidos.

Explorou-se um limite por inversa aproximada B: com `E=I-A*B` e `||E||F<1`,
`cond_2(A) <= ||A||F*||B||F/(1-||E||F)` em aritmética exata. Em N=320/1280,
o limite foi aproximadamente 806/2773 contra condições 163/327. Em N=1280,
inversa e residual levaram cerca de 0.27 s contra 0.73 s da SVD isolada, mas
exigem matrizes/produtos adicionais, são conservadores e não constituem limite
intervalar certificado em ponto flutuante. Em N=320 a alternativa foi mais
lenta. **Não foi adotada**: a revisão não introduz um novo gate aproximado ou
um backend LU validado, e mantém a SVD em toda avaliação.

## Testes finais, mudanças e pendências

| Execução | Aprovados | Falhas | Ignorados | Warnings pytest |
|---|---:|---:|---:|---:|
| Base v0.3.0, antes das mudanças | 337 | 0 | 21 | 0 |
| v0.3.1, comando padrão | 366 | 0 | 21 | 0 |
| v0.3.1, UIUC locais incluídos | 387 | 0 | 0 | 0 |

A base levou 65.11 s; as execuções finais levaram 66.51 s no modo padrão e
67.77 s com UIUC. Conferiram-se 51 artefatos, 43 hashes de entrada dos casos
e seis snapshots de código do lote final. A figura foi inspecionada. Um
`SolverResult` v0.3.0 previamente salvo também foi carregado com a v0.3.1.

São 28 testes novos em `test_panel_review.py` e um cenário adicional de parser
no arquivo existente: **29 casos novos**. Os 21 ignorados do modo padrão são
apenas as referências UIUC não redistribuídas. A suíte com dados locais verifica
seus hashes, preservando os contratos e testes geométricos anteriores.
Também se verifica a referência independente com produção indisponível,
quadratura, Blasius, refinamento, falha de condicionamento, equivalência de
montagem, cancelamento distante, flag obrigatória de fuga finita, serialização
e hashes/proveniência da CLI. A execução da CLI usa `AttemptContext.run_process`.

Arquivos de produção alterados: kernel de influência/montagem, identidade e
diagnóstico do adaptador, versões. Não se alteraram `AerodynamicSolver`,
`OperatingConditions`, runner, coeficientes, referência de momento, esquemas
neutros ou controle de processos. `aeroopt.panel-output.v1` ganha o campo de limitações
da versão do adaptador; resultados neutros antigos continuam carregáveis.
Documentação, contexto persistente e exemplos da revisão foram atualizados.

Pendências explícitas: comparação numérica com benchmark publicado exato
multielemento; precisão de Cp junto à cúspide; CM pequeno em limites práticos;
fechamento/cantos de fuga finita; resolução de slots estreitos e faixas de
geometria/incidência; custo de campanhas incluindo auditoria e refinamento.
Não se implementaram camada limite, transição, turbulência, arrasto viscoso,
estol, CFD geral ou algoritmos de otimização. **Software passando nos testes,
verificação numérica e validação física são níveis distintos; nenhum dos dois
primeiros estabelece aptidão para otimização de desempenho real.**
