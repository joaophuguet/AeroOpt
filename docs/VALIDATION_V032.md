# AeroOpt v0.3.2 — benchmark independente e confiabilidade numérica

## Estado científico

Esta etapa acrescenta uma solução analítica independente de dois corpos,
comparações quantitativas com Hess–Smith (HS) e vorticidade linear (LV), e
diagnósticos explícitos para ordenação numérica de candidatos. **A reconstrução
dos aerofólios publicados de Williams não foi certificada. O critério de
aceitação do benchmark publicado independente permanece não satisfeito.**
Os círculos de Williams foram reproduzidos na precisão impressa; o caso
conforme completamente declarado fornece evidência adicional, sem se fazer
passar pela geometria original. Não há validação experimental ou qualificação
para otimização aerodinâmica automática.

| Nível | Evidência e alcance |
|---|---|
| Correção matemática | Reflexões em círculos, circulação, derivadas conformes, sinais e normalização explicitados; resíduos de contorno/estagnação e duas integrações analíticas independentes conferidos numericamente. Não é uma prova formal do código. |
| Verificação numérica | Refinamento HS/LV contra o escoamento analítico da geometria declarada; erros de CL, CM, circulação e Cp; gates discretos mantidos. |
| Concordância independente publicada | Circulações e pontos de estagnação dos **círculos** concordam com a Tabela 1; a reprodução do aerofólio final falha na precisão publicada. Aceitação correspondente pendente. |
| Validação física | Não realizada. Todos os casos são potenciais, 2D, incompressíveis e invíscidos. |
| Otimização preliminar | Apenas evidência local de ordenação numérica; nenhuma qualificação para busca automática de desempenho real. |

## Fonte, proveniência e leitura

B. R. Williams, *An Exact Test Case for the Plane Potential Flow About Two
Adjacent Lifting Aerofoils*, ARC R&M 3717, relatório de 1971, publicação de
1973. [PDF original](https://reports.aerade.cranfield.ac.uk/bitstream/handle/1826.2/2993/arc-rm-3717.pdf?isAllowed=y&sequence=1).
O endereço canônico devolveu HTTP 403 ao download; o mesmo arquivo foi obtido
pelo caminho direto `bitstream/1826.2/2993/1/arc-rm-3717.pdf`, conservado
integralmente em `outputs/reference_data_v032/arc-rm-3717.pdf`:

```text
SHA-256 db8ae04e423bd17becf7a1615db477246a576d4dd8300bc6f1dcfdd3c5b9a369
1 349 341 bytes
```

`provenance.json` registra URL, data de obtenção, tamanho e hash. O relatório
foi extraído e suas 38 páginas renderizadas; conferiram-se visualmente as
equações de mapeamento, apêndices I–IV e tabelas selecionadas. Os utilitários
MiKTeX existentes falharam por instalação inicial incompleta, preservando logs.
A leitura usou PyMuPDF isolado em `.tmp/pdf_tools_v032`, sem dependência do
projeto ou instalação global. Texto e imagens ficam em
`outputs/reference_data_v032/inspection_pymupdf/`.

A transcrição manual seletiva `published_checks.json` conserva os valores
arredondados das Tabelas 1–3 e sua fonte; o catálogo versionado
`examples/reference_catalog_v032.json` fixa seu SHA-256
`b937b1e022e7eb789cba311465625062cb4294beb37a943f76c57a556bee9778`.
Não são coordenadas contínuas exatas nem uma digitalização ajustada. O PDF e
os pontos externos permanecem em `outputs/`, fora do Git, pois não foram
estabelecidas condições de redistribuição. O catálogo e parâmetros
matemáticos do problema são versionados.

O relatório descreve: §2, superposição dos escoamentos uniforme e circulantes
em dois círculos; §3, dois mapeamentos e regra da cadeia para velocidades;
§4 e Tabelas 1–3, configurações A/B e pressões/forças; Apêndice I, imagens de
dupleto/vórtice; II, pontos inversos complementares; III, convergência dos
dupletos; IV, agrupamento de vórtices opostos que torna convergente a série.
As figuras auxiliam a interpretação dos eixos; não forneceram pontos para o
solver. As forças por elemento não podem ser inferidas isoladamente de sua
circulação em um escoamento com interação. O relatório não fornece CM para
a referência de quarto de corda adotada aqui.

## Solução independente e convenções

`examples/_williams_reference.py` não importa nenhum dos kernels de painéis.
Usa velocidade complexa `w=u-iv`, velocidade livre unitária `exp(-i*alpha)`
e circulação **CCW positiva**, oposta à circulação positiva impressa nas
tabelas de Williams. Para um círculo de centro C e raio R:

```text
I_C(p) = C + R² / conj(p-C)
```

O dupleto original tem intensidade `D=R²*conj(U)` e contribuição de velocidade
`-D/(z-C)²`. A imagem em outro círculo tem
`D'=-R'²*conj(D)/conj(p-C')²`, na posição `I_C'(p)`.
As duas famílias de imagens são somadas com 64 reflexões alternadas. Um
vórtice unitário CCW contribui `1/(2*pi*i*(z-C))`; cada imagem acrescenta um
vórtice oposto no inverso e um positivo no centro. Os pares opostos são
avaliados como `(p-q)/(2*pi*i*(z-p)*(z-q))`, evitando subtrair termos grandes
que isoladamente não decrescem. Isso implementa o agrupamento físico do
Apêndice IV. Não se avalia um potencial logarítmico multivalorado de vórtices;
suas derivadas são unívocas no fluido. As circulações são obtidas resolvendo
um sistema 2×2 de velocidade tangencial nula nos pontos prescritos de fuga.

Os círculos têm raios `a=1.096`, `b=.5` e deslocamento `1.9-.4i`. Para o
primeiro mapeamento, `C1=1-a*exp(-i*.001)`, `C2=C1+1.9-.4i`, `T1=1` e
`T2=C2+b*exp(-i*(.08725+eta))`. A corrente livre é `alpha_raw=0` nesse
referencial de mapeamento. **As coordenadas de estagnação impressas na
Tabela 1 estão no referencial alinhado com a linha dos centros**: é necessário
rotacionar por `exp(-i*arg(1.9-.4i))`. Confundir esses eixos altera incidência
e circulações; não se usou essa rotação para alterar a geometria de entrada.

A transformação adotada é:

```text
K(z;c,n) = n*c*(1+r**n)/(1-r**n), r=(z-c)/(z+c)
n = 2 - .17279/pi
K'(z) = 4*n²*c²*r**n / ((z²-c²)*(1-r**n)²)
```

`r**n` usa `exp(n*Log(r))` principal. O corte `[-c,c]` deve ficar dentro do
corpo sólido. Os testes verificam sua contenção nos dois estágios e a
regularidade das poligonais. Nos pontos críticos ±c, usam-se os limites exatos
`K=±n*c`, `K'=0`. Os pontos de fuga são críticos e de estagnação;
`Cp_TE=1` é o limite para esta **cunha de ângulo não nulo e espessura zero**,
sem deslocamento epsilon. Isso não se aplica à cúspide de Joukowski ou à
tampa de espessura finita.

### Obstáculo de reconstrução e alternativa declarada

Não se conseguiu resolver a interpretação do centro/escala do segundo
mapeamento, após a primeira transformação que distorce o segundo círculo,
de forma a reproduzir as tabelas na precisão impressa. A descrição de §3,
`k2=1.15` e a rotação eta foram investigadas sem ajustar coeficientes para
forçar concordância. A limitação é da reconstrução implementada; não se
afirma que o problema publicado seja mal definido ou impossível de reproduzir.

A alternativa fixa explicitamente, com `m=K1(z)`:

```text
R2m = |K1(T2)-K1(C2)|
c2 = R2m / 1.15
O2 = K1(T2) - c2*exp(-i*eta)
F(z) = O2 + exp(-i*eta)*K2(exp(i*eta)*(K1(z)-O2);c2,n)
```

Essa escolha não é inferida como exata do relatório. Define uma geometria
própria, reproduzível, e um escoamento independente por `w_F=w_z/F'(z)`.
As duas configurações chamam-se **A_like/B_like**, nunca A/B certificadas.
O caso simples de dois círculos, anterior ao mapeamento, fornece o controle
direto independente que efetivamente coincide com dados publicados.

O LE nominal é a imagem da interseção esquerda do círculo principal com o
eixo real do primeiro mapa: `z_LE=C1-a*exp(+i*.001)`. Não é mínimo x
amostrado, ponto mais distante da fuga nem ponto diametralmente oposto a T1.
Com `L=F(z_LE)` e `d=F(T1)-L`, as entradas são `(F-L)/d`;
`alpha_normalized=alpha_raw-arg(d)`. LE e TE principais ficam exatamente
em `(0,0)` e `(1,0)`, com índice explícito `N/2` do LE. A comparação com a
tabela usa adicionalmente `(F-L)/real(d)`, conservando o pequeno y da fuga
impresso no relatório; essa convenção não é ajustada aos pontos externos.

Cada malha analítica é gerada explicitamente a partir do parâmetro angular
dos círculos, em sentido CCW a partir da fuga. No principal, os dois
intervalos encontram o LE nominal em `theta=pi+.002`; no flap a distribuição
é uniforme. Não se rediscretizam entradas externas. Ambos os métodos recebem
os mesmos arrays, incidência e referência de momento `(0.25,0)`, CCW,
normalizado por `q_inf*c_main²`. A referência integra Cp analítico por 8192 e
16384 pontos/corpo; Blasius num contorno exterior circular de raio 8 com 1024
pontos fornece uma integração independente de força/momento. A identidade
de força total usa `CL=-2*sum(Gamma_normalized)`. Erros de truncamento da
série, quadratura e discretização são registrados separadamente.

## Reprodução e artefatos

```powershell
.\.venv\Scripts\python.exe examples/verify_ranking_reliability.py --williams-data outputs/reference_data_v032 --output outputs/verification_v032_final
.\.venv\Scripts\python.exe -m pytest -q
.\.venv\Scripts\python.exe -m pytest -q --reference-data=outputs/reference_data_v011 --williams-data=outputs/reference_data_v032
```

O diretório de saída deve estar novo/vazio e dentro do workspace. O comando
também funciona offline sem `--williams-data`: mantém a referência declarada
e ignora apenas as comparações com o PDF/transcrição locais. `--quick` é um
smoke test, sem valor de campanha completa. A transcrição usada aqui não é
redistribuída; obter somente o PDF não fornece automaticamente esse sidecar.
Os hashes e critérios para uma transcrição equivalente estão no catálogo e
na função `load_published`.

O lote conserva inputs analíticos completos, hashes, Cp, circulações, forças,
momentos, diagnósticos discretos, referências e tempos em JSON; tabelas CSV e
figura; tentativas HS da campanha NACA executadas pelo runner; cópia local dos
originais; snapshots do código de pesquisa e kernel; inventário SHA-256 de
todos os artefatos. O inventário não inclui seu próprio hash. Os arquivos
derivados de dados externos permanecem em `outputs/`.

## Ordenação numérica e custos

As regras completas e API estão em [NUMERICAL_RELIABILITY.md](NUMERICAL_RELIABILITY.md).
O indicador por candidato é o maior entre o último incremento de malha de
cada método e a discrepância entre métodos na malha mais fina. A margem de
comparação é `2*(U_A+U_B)`. Exigem-se refinamento contrativo, não oscilatório,
ausência de limitações e sinal comum das diferenças nos dois níveis mais
finos de ambos os métodos. Uma troca de sinal é inconsistente; separação
insuficiente ou evidência ausente é inconclusiva. Não se ajustou a política
após observar os resultados, nem se extrapolou um limite universal de erro.

Campanha fixa: NACA 0012 principal, NACA 4412 flap, fugas fechadas, razão de
cordas .3, deflexão 25°, LE do flap `(1.02,-.04)`, alpha 6°. Variações
individuais: deflexão ±1°, x ±.01, y ±.005, alpha ±.5°; perturbações próximas
de +.001°, +.00001 em x/y e +.0001° em alpha. São 13 candidatos e 78 pares
por coeficiente; **nenhuma busca ou otimização**. As quatro malhas têm
80/160/320/640 painéis por elemento (160 a 1280 totais). O runner exige
folga .002; as geometrias passam pela validação a cada nível. Reynolds
200000 é metadado explícito do runner, sem efeito no modelo invíscido;
Mach=0. A referência de momento é explícita. LV permanece um código de
pesquisa, com seus próprios gates de condição e resíduos.

Os tempos HS incluem runner, auditoria, serialização e disco; os tempos LV
medem apenas a chamada ao solver de pesquisa. Não são uma comparação pura
de kernels. O tempo de referência analítica é separado. Todos os custos são
medições locais, dependem de BLAS/hardware/I/O, não são promessas de desempenho.

## Investigações adicionais e decisões

A fuga finita foi estudada mantendo a mesma poligonal e subdividindo apenas
os segmentos retos da tampa em 1/4/8/16 partes, nas quatro malhas de superfície.
Assim se distingue subdivisão local de mudança de geometria. A condição de
Kutta continua aplicada às superfícies terminais superior/inferior; resolver
essa equação discretizada não prova validade do modelo nos cantos da tampa.
O benchmark de fuga afiada com cunha não valida uma tampa finita. Não se
introduziram pressão de base empírica, viscosidade ou correção de força.
`FINITE_TE_UNVERIFIED` permanece obrigatório e impede qualificação de rankings.

Joukowski foi investigado em 160/640/1920 painéis, com a referência de momento
analítica da v0.3.1. `examples/_pressure_quadrature.py` calcula integração
Gauss de ordem 2/4 no campo HS já resolvido. Explicita a auto-influência no
interior de cada painel e registra a velocidade normal fora dos pontos de
colocação. Avalia Cp tangencial e Cp de velocidade total como diagnósticos.
Não resolve impermeabilidade nesses novos pontos. Reduzir erro integrado
assim não basta para justificar substituir a formulação: os resíduos normais
fora da colocação e a sensibilidade à ordem de quadratura são significativos.
Nenhuma dessas saídas é publicada como coeficiente de produção; nenhuma
correção foi adotada. A limitação de CM pequeno na cúspide permanece.

## Compatibilidade

Não se alteraram influências HS, Kutta, gates, `AerodynamicSolver`, runner,
schemas neutros, sinal de momento, coeficientes ou controle de processos.
O pacote passa a 0.3.2; **a identidade científica do adaptador permanece
0.3.1**, pois a implementação de produção não mudou. Não se invalida o hash
de análises idênticas por acrescentar ferramentas de verificação externas.
`cd=None` e os avisos existentes permanecem. O novo pacote
`aeroopt.verification` é optativo e separado dos coeficientes físicos.

## Resultados quantitativos do lote final

Lote local `outputs/verification_v032_final/`, Python 3.13.15, NumPy 2.5.4,
Windows, ramo `feature/panel-solver-v0.3`, base v0.3.1
`e2eb380cf9796f0612eb020d503be708d5621cd2`. Todas as soluções HS/LV da
campanha foram aceitas pelos respectivos gates. O lote levou 125.07 s até
concluir os cálculos; a execução supervisionada, incluindo importação,
figura, snapshots e inventário, levou 130.49 s.

### Controles analíticos e discrepância publicada

| Grandeza | A_like | B_like |
|---|---:|---:|
| Circulações dos círculos, CCW | -5.169208727 / -1.778027647 | -3.153751732 / -.655201646 |
| Máx. diferença frente à circulação impressa dos círculos | 2.0873e-4 | 2.4827e-4 |
| Máx. diferença de coordenada dos pontos de estagnação no eixo dos centros | 4.7975e-5 | 3.0019e-5 |
| Máx. velocidade normal nos círculos | 2.881e-15 | 1.704e-15 |
| Máx. velocidade nos pontos prescritos de estagnação | 3.001e-16 | 6.825e-16 |
| Diferença de Gamma com 32/64 reflexões, em ponto flutuante | 0 | 0 |
| Incidência no referencial normalizado, graus | -.3376142127 | 1.5003037663 |
| CL por circulação, referência declarada | 3.737898256341 | 2.028017349646 |
| CM por Blasius, CCW em .25c | 1.265482463030 | .620166229855 |
| Gamma/c principal / flap | -1.390623209481 / -.478325918690 | -.839582766405 / -.174425908418 |
| Diferença de CL por pressão, 8192 → 16384 pontos/corpo | 1.082e-9 | 1.307e-9 |
| Diferença de CM por pressão, mesmo refinamento | 9.445e-10 | 9.917e-10 |
| CL por pressão a 16384 menos CL por circulação | 4.090e-10 | 4.769e-10 |
| CM por pressão a 16384 menos CM por Blasius | 3.549e-10 | 3.618e-10 |

As circulações dos círculos publicadas, positivas CW, são `5.169/1.778` para
A e `3.154/.655` para B. Os limites de arredondamento são .0005 para essas
circulações e .00005 para coordenadas de quatro casas; os controles passam
sem afrouxar a precisão impressa. Isso não constitui reprodução do aerofólio.

No aerofólio transformado, o máximo erro de circulação normalizada frente à
tabela é `2.768e-4` em A_like e `4.172e-4` em B_like, acima de `.00005`.
As diferenças de CL total por circulação são `-.00070174` e `-.00098265`.
Nos dez pontos selecionados de A, a maior distância do ponto impresso à curva
declarada é `.00028125`, maior que o arredondamento de coordenadas de cinco
casas. Na fuga do flap, por exemplo, a distância é `.00027362`.
Essas discrepâncias foram **conservadas**. A projeção
mais próxima é apenas diagnóstico dos valores arredondados e não foi usada
para gerar, reparar ou ajustar entradas. Os Cp projetados e suas diferenças
estão no JSON, sem serem apresentados como erros de discretização contra o
aerofólio publicado, pois as geometrias não coincidem. Não há comparação de
CM publicado, pois esse valor não é fornecido para a referência adotada.

### Erros frente ao escoamento analítico da geometria declarada

N é o número **total** de painéis, dividido igualmente entre os elementos.
Os erros de CL/CM são assinados, numérico menos referência; Gamma mostra o
máximo módulo de erro por elemento. Cp é RMS por elemento, principal/flap.
CM usa Blasius como referência. Nenhum valor é uma correção aplicada ao solver.

| Caso | N | Método | Erro CL | Erro CM | Máx. erro Gamma | RMS Cp principal/flap |
|---|---:|---|---:|---:|---:|---|
| A_like | 160 | HS | -.0652004 | -.0451419 | .0134950 | .131802 / .082162 |
| A_like | 320 | HS | -.0308170 | -.0225803 | .00533826 | .068640 / .038837 |
| A_like | 640 | HS | -.0144136 | -.0110966 | .00206167 | .035413 / .018045 |
| A_like | 1280 | HS | -.00675671 | -.00542830 | .000784536 | .017954 / .008461 |
| A_like | 160 | LV | -.0348923 | -.00729181 | .00230137 | .167261 / .155000 |
| A_like | 320 | LV | -.0167764 | -.00360579 | .000637486 | .082801 / .076109 |
| A_like | 640 | LV | -.00817739 | -.00177987 | .000170040 | .042759 / .037311 |
| A_like | 1280 | LV | -.00402789 | -.000881745 | .0000437668 | .022868 / .018772 |
| B_like | 160 | HS | -.0367701 | -.0264562 | .00798580 | .053203 / .033443 |
| B_like | 320 | HS | -.0168662 | -.0130202 | .00292910 | .027597 / .016592 |
| B_like | 640 | HS | -.00775946 | -.00637028 | .00103776 | .014182 / .008278 |
| B_like | 1280 | HS | -.00360133 | -.00311700 | .000351765 | .007178 / .004149 |
| B_like | 160 | LV | -.0207355 | -.00418173 | .00187194 | .066644 / .021402 |
| B_like | 320 | LV | -.00977616 | -.00200516 | .000490813 | .033462 / .012322 |
| B_like | 640 | LV | -.00472670 | -.000979229 | .000125471 | .017670 / .007682 |
| B_like | 1280 | LV | -.00232122 | -.000483621 | .0000313841 | .009727 / .004983 |

Os módulos de todos esses erros decrescem nas quatro malhas. Para CL/CM, a
razão de erros entre as duas malhas mais finas é aproximadamente 2, coerente
com convergência próxima de primeira ordem nesses casos. Não se extrapolou
essa taxa para outros perfis ou coeficientes. Em A_like, HS chega a cerca de
.181% de erro CL e .429% de erro CM em N=1280; LV tem menor erro integrado e
de circulação, porém maior RMS de Cp. Concordância integrada não garante Cp
local de maior qualidade, particularmente junto à fuga e regiões de sucção.

O Cp exato é amostrado no ponto médio do parâmetro da curva; o Cp numérico
está no centro do segmento reto. A diferença entre locais, parte do erro de
geometria/discretização, fica explícita: em N=1280 seus máximos são
`6.850e-6/4.219e-6` em A_like e `9.089e-6/4.099e-6` em B_like. Não se
interpolou a solução numérica ou a tabela externa. Não há diferença de modelo
físico entre os dois solvers e a referência declarada: todos usam o mesmo
problema potencial. Há diferenças de representação do contorno, singularidades,
Kutta discreto e quadratura de pressão.

Os máximos `cond_2` do benchmark são 457.72 (HS) e 708705.86 (LV), ambos abaixo
de 1e10. Os resíduos máximos de Kutta são `3.331e-15` (HS) e zero em ponto
flutuante (LV, soma algébrica das vorticidades terminais); os normais são
`6.323e-14` e `1.743e-14`. Não são a mesma equação de Kutta: HS equilibra
velocidades das superfícies terminais e LV impõe vorticidades nodais opostas.
Os resíduos pequenos coexistem com erros de CL/CM não desprezíveis. Não se
usa resíduo como estimativa do erro do contínuo. Todos os valores e painéis
terminais estão nos JSON dos casos. A figura `verification.png`, inspecionada,
mostra erros integrados de A_like e as contagens de classificação; B_like e
Cp/Gamma são conservados na tabela CSV e dados completos.

### Resultados de ordenação

| Coeficiente | Pares distinguíveis | Inconclusivos | Inconsistentes | Total |
|---|---:|---:|---:|---:|
| CL | 54 | 24 | 0 | 78 |
| CM | 42 | 36 | 0 | 78 |

Não apareceram trocas de sinal nos **dois níveis mais finos** deste conjunto.
Isso não prova estabilidade em todos os níveis ou fora desses candidatos.
Os testes incluem deliberadamente trocas de ordenação entre métodos/malhas
para verificar que são classificadas como inconsistentes. Todas as histórias
da campanha têm últimos incrementos não oscilatórios e decrescentes em módulo.

História da base:

| N total | CL HS | CL LV | CM HS | CM LV |
|---:|---:|---:|---:|---:|
| 160 | 3.14061203 | 3.11994344 | .79406030 | .80661286 |
| 320 | 3.14942923 | 3.13855703 | .80195983 | .80810754 |
| 640 | 3.15316882 | 3.14748734 | .80591799 | .80892273 |
| 1280 | 3.15478454 | 3.15185772 | .80787128 | .80934943 |

Na base, U_CL=.00437038 e U_CM=.00195329. Não são limites do erro real.
Os exemplos abaixo mostram intervalos de diferenças B−base nos dois níveis
mais finos dos dois métodos e a margem heurística adotada:

| B | Coef. | Intervalo B−base | Margem | Resultado |
|---|---|---|---:|---|
| deflexão +1° | CL | [.0671330,.0673625] | .0177291 | distinguível, B maior |
| x +.01 | CL | [-.0164634,-.0161801] | .0175570 | inconclusivo |
| y +.005 | CL | [.0209055,.0210921] | .0174569 | distinguível, B maior |
| alpha +.5° | CL | [.0696053,.0698998] | .0178057 | distinguível, B maior |
| deflexão +.001° | CL | [6.7484e-5,6.7711e-5] | .0174818 | inconclusivo |
| x +.00001 | CL | [-1.3709e-5,-1.3363e-5] | .0174816 | inconclusivo |
| y +.00001 | CL | [4.2597e-5,4.2944e-5] | .0174815 | inconclusivo |
| alpha +.0001° | CL | [1.3946e-5,1.4004e-5] | .0174816 | inconclusivo |
| deflexão +1° | CM | [.0217193,.0217552] | .00788476 | distinguível, B maior |
| x +.01 | CM | [.0128663,.0131477] | .00764446 | distinguível, B maior |
| y +.005 | CM | [.0008661,.0008956] | .00775994 | inconclusivo |
| alpha +.5° | CM | [.0028551,.0029146] | .00784681 | inconclusivo |

As quatro perturbações mínimas são inconclusivas também em CM. Mesmo quando
ambos os métodos concordam no sinal, a política não certifica diferenças muito
menores que os indicadores. Um aumento algébrico de CM não é uma melhoria de
estabilidade ou desempenho. A comparação de candidatos não pode herdar
automaticamente a precisão observada no benchmark A_like/B_like.

### Tempos

Medianas de 13 configurações por nível; o HS inclui runner/I/O, LV apenas
cálculo. As histórias incluem os dois métodos nas quatro malhas.

| N total | HS com runner, s | LV, s |
|---:|---:|---:|
| 160 | .416 | .061 |
| 320 | .520 | .231 |
| 640 | .722 | .749 |
| 1280 | 1.624 | 3.016 |

Uma história de quatro malhas custa mediana **7.488 s/configuração**, intervalo
7.067–7.815 s. A referência analítica com duas quadraturas e Blasius levou
.083 s em A_like e .157 s em B_like, fora dos solvers. Um resultado HS de
1280 painéis custa cerca de 1.62 s nessa campanha; comprar evidência de
refinamento e comparação independente custa mais. Não existe resolução
universalmente suficiente. Mesmo com esse custo, os casos mínimos continuam
inconclusivos. Não se aciona essa campanha a cada chamada de produção.

### Fuga finita e momento de Joukowski

Valores de CL para a mesma geometria finita, variando só as subdivisões da
tampa; N_s é a quantidade de painéis de superfície por elemento, sem a tampa:

| N_s | Tampa 1 | Tampa 4 | Tampa 8 | Tampa 16 | CM, tampa 1 → 16 |
|---:|---:|---:|---:|---:|---|
| 80 | 3.142203 | 3.164179 | 3.165972 | 3.166668 | .803054 → .810045 |
| 160 | 3.142743 | 3.167403 | 3.170059 | 3.170979 | .804701 → .812643 |
| 320 | 3.131602 | 3.161500 | 3.164805 | 3.166248 | .802258 → .811967 |
| 640 | 3.118663 | 3.158319 | 3.162336 | 3.164007 | .798983 → .811717 |

No nível mais fino, trocar 1 por 16 partes muda CL em .0453442 e CM em
.0127339, embora o máximo resíduo de Kutta nesses ensaios seja `2.887e-14`.
Mesmo 8 → 16 ainda muda CL em .0016713. A sequência de refinamento da
superfície não estabelece um limite comum com precisão certificada. Subdivisão
local não resolveu o modelo de fuga finita: flag e restrições foram mantidas.

Joukowski a alpha=4°, CM analítico CCW de quarto de corda
`.00124065623251238`:

| Painéis | CM por ponto médio (produção) | CM Gauss 2, pesquisa | CM Gauss 4, pesquisa | Máx. velocidade normal Gauss 2 / 4 |
|---:|---:|---:|---:|---|
| 160 | -.0152766074 | -.0083528746 | -.0046904090 | .07436 / .25364 |
| 640 | -.0050588218 | -.0021283592 | -.0008151304 | .07938 / .28353 |
| 1920 | -.0014082019 | -.0001140986 | .0004209244 | .08058 / .29101 |

Em 1920, o erro absoluto de CM cai de .00264886 para .000819732 no diagnóstico
Gauss 4, mas ainda é aproximadamente 66% do pequeno momento exato. Os resíduos
normais fora da colocação não são próximos de 1e-9 e não decrescem com essa
quadratura; aumentar a ordem muda significativamente o resultado. O Cp de
velocidade total dá valores próximos, sem resolver a inconsistência de
contorno. Isso apoia investigar representação/colocação/Kutta de ordem maior,
não publicar uma correção automática a partir da melhora integrada observada.
Produção e limitações foram preservadas.

## Regressão, auditoria e pendências

| Execução | Aprovados | Ignorados | Falhas | Tempo, s |
|---|---:|---:|---:|---:|
| Base v0.3.1, antes das alterações | 366 | 21 | 0 | 76.75 |
| v0.3.2, modo padrão | 403 | 23 | 0 | 105.61 |
| v0.3.2, UIUC e Williams locais | 426 | 0 | 0 | 85.61 |

Não houve warnings pytest. São 39 casos novos, incluindo dois opcionais com
os originais Williams: os 23 ignorados no modo padrão são 21 UIUC e dois
Williams. O modo padrão também executa integralmente a referência analítica
declarada e sua CLI, sem internet ou dados externos. O maior tempo da primeira
execução não estabelece uma regressão do kernel: a campanha da CLI acrescenta
trabalho e os tempos locais variam. O pacote foi reinstalado editável em .venv
usando apenas o wheel de build já disponível no workspace; `pip check` passou,
sem nova dependência de execução. Não se recriou a .venv nem houve instalação
global. O resumo de regressão foi preservado junto aos logs de lançamento.

Os novos testes cobrem fonte local fixada por hash, referência independente
com kernels indisponíveis, eixos/LE/TE/incidência, cortes/derivadas/limites,
covariância por rotação, truncamento das séries, Blasius/pressão/circulação,
refinamento de CL/CM/Gamma/Cp, Kutta, ausência/falha/limitações, ordenações
resolvidas/inconclusivas/inconsistentes e política; a CLI preserva hashes,
inputs, inventário e proteção contra sobrescrita. O subprocesso da CLI usa
`AttemptContext.run_process` com timeout e controle de descendentes.

Foram conferidos **646 artefatos** por SHA-256 e os oito snapshots contra as
fontes usadas no lote; a figura foi inspecionada. Permanecem pendentes a
reconstrução certificada do aerofólio publicado, resolução geral de slots e
geometrias/incidências diferentes, Cp junto a singularidades, CM pequeno de
cúspide e o fechamento finito. A validade física exige evidência adicional
compatível com Reynolds, viscosidade, transição, separação e estol, além da
conexão ao problema 3D de desempenho da aeronave. **Recomendação: manter o
uso em estudos invíscidos controlados; não iniciar otimização aerodinâmica
automática como se a precisão e a física estivessem qualificadas.**

Hash do inventário final `artifacts.json`:
`19797d7aaf1dbb7cd3f9b2f9996ebbfab9f27a01d30166779aaa7f53bf079435`.
