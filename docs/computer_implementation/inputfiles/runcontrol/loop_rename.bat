setlocal ENABLEDELAYEDEXPANSION  
set list=card1 card10 card11 card11a card11b card12 card12a card13 card14 card14a card15 card16 card17 card18 card19 card2 card20 card21 card22 card22b card23 card24 card25 card26 card27 card28 card29 card3 card30 card31 card32 card33 card34 card35 card36 card36a card36b card37 card38 card39 card4 card40 card41 card42 card42a card43a card43b card43c card43d card43e card44 card45 card45a card45b card45c card45d card46 card46a card46c card46e card47 card48 card49 card5 card50 card51 card52 card53 card54 card55 card56 card57 card58 card59 card6 card60 card61 card62 card63 card64 card65 card66 card66a card66b card67 card68 card69 card7 card70 card71 card71a card71b card72 card73 card74 card75 card76 card77 card78 card79 card8 card80 card81 card82 card83 card84 card85 card86 card87 card88 card89 card9 card90 card91 card91a card91b card92a
set ext = .rst
for  %%a in (%list%) do (
	@echo .. %%a: >>%%a.rst
	@echo. >> %%a.rst
	@echo %%a >>%%a.rst
	@echo ------- >> %%a.rst
	@echo. >> %%a.rst
	@echo  .. literalinclude:: ../efdc_card_opts/%%a  >> %%a.rst
)